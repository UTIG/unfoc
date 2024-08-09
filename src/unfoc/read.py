#!/usr/bin/env python3


"""
read.py

This module provides functions and objects for reading UTIG
radar data formats into numpy.

Possible input data formats include RADnh3/RADnh5 breakout data (e.g., HiCARS2, MARFA),
as well as HiCARS1 (RADjh1).

For RADnh3/RADnh5, we assume that the standard breakout procedure has been run
on the ELSA data files to produce our standard bxds and ct files.

NB: In theory we could write a parser to directly parse ELSA data files and provide
the same data ingest capability.

For RADjh1 data, we assume that the standard breakout procedure has been run on
the HiCARS1 data files and that the standard bxds1, bxds2, and ct files are available
in one directory.


The easiest way to read raw radar data is to use the RadBxds and RADjh1Bxds classes.

The most memory-efficient way to read RADnh3/RADnh5 data is to use the read_RADnhx_gen
generator.  This omits the step of pre-walking the data file to find radar records.

If you don't want to read all the radar data, you can use the index_RADnhx_gen function
to get an index of radar records present.

Right now, the RADjh1 and RADnh3/RADnh5 functions are a bit disjointed but in the
future you can expect them to be a little more unified and for the module to detect
the format of the data from the binary contents of the file.



"""

from collections import namedtuple, defaultdict
import itertools
import logging
import os
import struct
import gzip
import mmap

import numpy as np

# Pik1ChannelSpec
#import unfoc.parse_channels

__version__ = '2.1.0'

HEADER_FORMATS = {
    'RADnh3': {
        'header_t': namedtuple('radnh3_header',
                              'nsamp nchan vr0 vr1 choff ver resvd2'),
        'fmtstr':  '>HBBBBBB',
    },
    'RADnh5': {
        'header_t': namedtuple('radnh5_header',
                              'nsamp nchan vr0 vr1 choff ver resvd2 absix relix xinc rseq scount tscount'),
        'fmtstr': '>HBBBBBBddfLHL',
    }
}



class Trace:
    """
    Pairs channel + data and CT metadata for a trace.
    Data is either raw amplitudes, or complex amplitudes
    """
    def __init__(self, channel, data, ct):
        # type: (int, np.ndarray, Union[tuple,List[tuple]]) -> None
        self.channel = channel
        self.data = data # type: np.ndarray
        self.ct = ct # type: tuple




def get_radar_stream(filename:str):
    # type: (str) -> str
    # both RADnh3 and RADnh5 are the same first elements in header format
    # TODO: add support for detecting RADjh1
    fmtstr = '>HBBBBBB'
    fmtlen = struct.calcsize(fmtstr)
    with open(filename, 'rb') as fd:
        # read header
        buff = fd.read(fmtlen)
        if len(buff) < fmtlen: # pragma: no cover
            raise IOError("Unable to read %d bytes from %s." % (fmtlen, filename))

        # TODO: This check should only  happen once, not every read!!
        # Check it then set a flag and move the file pointer back to the start!
        # Check the version byte and see if it is RADnh3 or RADnh5
        try:
            v = ord(buff[6]) # python 2
        except TypeError:
            v = int(buff[6]) # python 3
        if v == 5:
            return "RADnh5"
        elif v == 0 or v == 255:
            return "RADnh3"
        else:
            raise ValueError("Can't determine stream for file %s" % filename)


CT_t = namedtuple('CT', 'seq tim')
def gen_ct(bxdsfile):
    """ Generate ct data from the ct file associated with bxdsfile
    We only return relevant information. We actually only want the seq and ct tim.

    """
    datadir = os.path.dirname(bxdsfile)
    ctfile1 = os.path.join(datadir, 'ct')
    ctfile2 = os.path.join(datadir, 'ct.gz')
    if os.path.exists(ctfile1):
        fh = open(ctfile1, 'rt')
    elif os.path.exists(ctfile2):
        fh = gzip.open(ctfile2, 'rt')
    else: #pragma: no cover
        return
    for line in fh:
        fields = line.rstrip().split()
        # e.g., THW PBA0a X66a 9644392 2020 02 02 03 22 59 70 2762000089
        # yield ct, seq
        yield CT_t(int(fields[3]), int(fields[11]))
    fh.close()



def get_radar_type(bxdsfile, nrecords=2000, stream=None):
    """ Inspect a raw datafile and detect what type of radar it came from
    returns 'HiCARS2' if it is a 1-antenna radar, or 'MARFA' if it is a 2-antenna
    radar.

    Inspect a raw datafile and detect all of the channels in this file
    """
    wanted_channels = [0, 2, 4, 64, 64+2, 64+4]
    MAX_CHANNELS = len(wanted_channels)

    choffs = set()
    gen = index_RADnhx_bxds_mmap_(bxdsfile, stream=stream)
    # fpos, headerlen, header.choff, header.nsamp
    for ii, (_, _, choff, _) in enumerate(gen):
        if choff == 0xff:
            choff = 0x00
        choffs.add(choff)
        if ii >= nrecords or len(choffs) >= MAX_CHANNELS:
            break
    gen.close()

    chans = []
    for choff in choffs:
        chans.extend([choff+1, choff+2])

    logging.debug("get_radar_type found channels: %s", sorted(chans))
    if len(choffs) >= 3:
        return 'MPOL', chans
    if len(choffs) >= 2:
        return 'MARFA', chans
    else:
        assert len(choffs) == 1
        return 'HiCARS2', chans

def sync_radar_start(bxdsfile: str, nrecords:int=2000, stream=None):
    """ Inspect a bxds file and figure out how many channels it has,
    and what the first complete rseq number is.  This allows us to
    make sure that all of the radar samples are synchronized between
    channels 

    Return value:
        - The record number of this record in the bxds file
        - The file position in the bxds file
        - The rseq number
    """
    # mapping of rseq number and the file position of that record
    rseq_to_fpos = defaultdict(list)

    gen1 = index_RADnhx_bxds(bxdsfile, stream=stream, full_header=True)
    gen2 = itertools.islice(gen1, 0, nrecords)
    for ii, (fpos, _, header) in enumerate(gen2):
        rseq_to_fpos[header.rseq].append((ii, fpos))
        if len(rseq_to_fpos) > 10:
            break

    # Number of channels is the max number of times we see an rseq value
    nchan = max([len(fposlist) for fposlist in rseq_to_fpos.values()])
    for rseq, fposlist in sorted(rseq_to_fpos.items()):
        if len(fposlist) == nchan:
            # Return the first rseq and first position where this was seen.
            return *fposlist[0], rseq

    raise ValueError("No radar records in %s" % bxdsfile)

def index_RADnhx_bxds(input_filename, stream=None, full_header=False):
    # type: (str) -> Generator[tuple]
    """ Read the positions of packets within a RADnh3 and RADnh5 bxds file
    and return these as a generator
    Routines can then use this to seek to the correct location in a file.
    # TODO: rework this function to always return the full header and make
    all callers grab the member methods
    """

    if stream is None:
        stream = get_radar_stream(input_filename)

    if stream not in HEADER_FORMATS: # pragma: no cover
        raise ValueError("Invalid stream type for file: %s" % input_filename)

    header_t = HEADER_FORMATS[stream]['header_t']
    fmtstr = HEADER_FORMATS[stream]['fmtstr']
    fmtlen = struct.calcsize(fmtstr)
    # In theory, it's faster to do this here and only compile the format string once.
    header_struct = struct.Struct(fmtstr)

    with open(input_filename, 'rb') as fd:
        fpos = 0
        while True:
            # read header
            fpos = fd.tell()
            buff = fd.read(fmtlen)
            if len(buff) < fmtlen:
                return
            #ctdata = next(genct) if genct else (None, None)
            header = header_t._make(header_struct.unpack(buff))
            headerlen = fmtlen
            if stream == "RADnh5":
                # We don't do anything with the timestamps yet, but we need
                # to skip over 'em to get to the radar data.
                if header.tscount > 0:
                    # timestamp count should be on the order of the number of stacks.
                    assert(header.tscount < 100)
                    #tstamps = struct.unpack_from('>{:d}d'.format(tscount), fd.read(8*tscount))
                    fd.seek(8*header.tscount, os.SEEK_CUR)
                    headerlen += 8*header.tscount

            input_samples = header.nsamp
            assert 0 < input_samples < 10000

            if full_header: # want to make this the default behavior at some point
                yield fpos, headerlen, header
            else:
                yield fpos, headerlen, header.choff, header.nsamp #, ctdata)
            # Skip the rest of this record without yielding any values
            fd.seek(4*input_samples, os.SEEK_CUR)



def index_RADnhx_bxds_mmap_(input_filename, stream=None, full_header:bool=False):
    # type: (str) -> Generator[tuple]
    """ Read the positions of packets within a RADnh3 and RADnh5 bxds file
    and return these as a generator
    Routines can then use this to seek to the correct location in a file.

    This version uses a slightly different syntax but it's functionally identical
    to index_RADnhx_bxds
    """

    if stream is None:
        stream = get_radar_stream(input_filename)

    if stream not in HEADER_FORMATS: # pragma: no cover
        raise ValueError("Invalid stream type for file: %s" % input_filename)

    header_t = HEADER_FORMATS[stream]['header_t']
    fmtstr = HEADER_FORMATS[stream]['fmtstr']

    # In theory, it's faster to do this here and only compile the format string once.
    header_struct = struct.Struct(fmtstr)
    fmtlen = header_struct.size
    with open(input_filename, 'rb') as fd:
        mmbxds = mmap.mmap(fd.fileno(), 0, access=mmap.ACCESS_READ)
        filelen = mmbxds.size()
        fpos = 0
        while fpos + fmtlen < filelen:
            # read header
            headerdata = header_struct.unpack_from(mmbxds, fpos)
            fpos_next = fpos + header_struct.size
            #ctdata = next(genct) if genct else (None, None)
            header = header_t._make(headerdata)
            headerlen = header_struct.size
            if stream == "RADnh5":
                # We don't do anything with the timestamps yet, but we need
                # to skip over 'em to get to the radar data.
                if header.tscount > 0:
                    # timestamp count should be on the order of the number of stacks.
                    assert header.tscount < 100
                    fpos_next += 8*header.tscount
                    headerlen += 8*header.tscount

            assert 0 < header.nsamp < 10000

            if full_header: # want to make this the default behavior at some point
                yield fpos, headerlen, header
            else:
                yield fpos, headerlen, header.choff, header.nsamp #, ctdata)
            # Skip the rest of this record without yielding any values
            fpos_next += 4*header.nsamp
            fpos = fpos_next
        mmbxds.close()



def read_RADnhx_gen(bxds_filename, channel, stream=None):
    """ Return a sequence of traces from only one channel, from a bxds
    channel offset is a one-based index.

    Takes a bxds filename as an argument.
    Reads from this file and a ct file alongside it.

    """
    assert 1 <= channel # one-based channel number
    # choff1 is used to calculate the byte offset of radar data within the packet.
    choff1 = (channel - 1) % 2
    # Filter for listed channel offset in radar record
    chfilter = (channel - 1) - choff1

    with open(bxds_filename, 'rb') as fd: # for reading traces
        ctgen = gen_ct(bxds_filename)
        radgen = index_RADnhx_bxds(bxds_filename, stream)
        for radinfo, ctinfo in zip(radgen, ctgen):
            fpos, headerlen, rchoff, nsamp = radinfo
            if rchoff == 0xff:
                rchoff = 0

            # 0x40 is the flag for along-track transmit polarization
            # 0x0f is the channel offset field
            # We want to ignore other bits.
            rchoff &= 0x4f

            if chfilter != rchoff:
                continue

            i0 = fpos + headerlen + choff1 * nsamp * 2
            fd.seek(i0, os.SEEK_SET)
            trace1 = np.fromfile(fd, dtype='>i2', count=nsamp) 

            yield Trace(channel, trace1, ctinfo)


class RadBxds:
    """ Reader for random access to traces in a RADnh3 or RADnh5 bxds as numpy arrays
    Future: add support for also auto-detecting RADjh1 bxds
    alternative name ideas:
    
    """
    def __init__(self, filename=None, channel=None, stream=None):
        """ Initialize the reader to access one channel's records """
        self.index_ = []
        self.mmbxds_ = None
        self.fd_ = None
        self.cts_ = None
        if filename is not None:
            self.open(filename, channel, stream)

    def __del__(self):
        self.close()

    def close(self):
        if self.mmbxds_ is not None:
            self.mmbxds_.close()
            self.mmbxds_ = None
        if self.fd_ is not None:
            self.fd_.close()
            self.fd_ = None
        self.cts_ = None
        self.index_ = []


    def open(self, filename, channel, stream=None):
        """ Open the bxds and load record index for bxds
        This ignores sequence numbers and assumes that there are no
        out-of-order radar records within a channel.
        In theory we could sort the records by sequence number.

        filename: bxds filename
        channel: one-based channel number
        stream: (optional) hint for stream type

        """

        assert channel >= 1
        self.channel0_ = channel - 1 # zero-based channel index
        self.trace_byteoffset_ = 2 * (self.channel0_ & 1)
        self.fd_ = open(filename, 'rb')
        self.mmbxds_ = mmap.mmap(self.fd_.fileno(), 0, access=mmap.ACCESS_READ)

        self.index_ = []
        filesize = self.mmbxds_.size()
        # Channel offset that we are expecting to filter for.
        choff = self.channel0_ - (self.channel0_ & 1)
        # Length of a record's trace data in bytes per sample
        bytes_per_samp = ((self.channel0_ & 1) + 1) * 2

        for ii, item in enumerate(index_RADnhx_bxds_mmap_(filename, stream=stream)):
            if item[2] == 0xff: # Replace choff if it is 0xff
                item = (item[0], item[1], 0x00, item[3])
            if item[2] == choff:
                lastbyte = item[0] + item[1] + bytes_per_samp * item[3]
                if lastbyte <= filesize:
                    self.index_.append(item + (ii,))

    def size(self):
        """ Return the number of traces for this channel """
        return len(self.index_)

    def __len__(self):
        """ Return the number of traces for this channel """
        return len(self.index_)


    def __getitem__(self, idx):
        """ Return data for selected radar records as an ndarray.

        Records to return can be selected using slice notation
        Example:
        bxdsmm = RadBxds(bxdsfile, channel=2)
        # Get the 6th trace for channel 2
        trace = bxdsmm[5]
        # Return a 5 x 3200 array for channel 2.
        traces = bxdsmm[5:10]
        # Return a portion of the array
        traces = bxdsmm[5:10, 200:400]
        """

        if isinstance(idx, slice): # 1D slice case
            indices = idx.indices(len(self.index_))
            ntraces = len(range(*indices))
            data = np.empty((ntraces, self.index_[0][3]), dtype=">i2")

            for ii, idxinfo in enumerate(self.index_[idx]):
                fpos, headerlen, _, nsamp, _ = idxinfo
                i0 = fpos + headerlen + self.trace_byteoffset_ * nsamp
                data[ii, :] = np.frombuffer(self.mmbxds_, dtype=">i2", offset=i0, count=nsamp)
        elif isinstance(idx, tuple):
            # multiple indices -- pass others down to ndarray
            data = self.__getitem__(idx[0])
            if isinstance(idx[0], slice):
                return data[:, idx[1]]
            else: # single integer
                return data[idx[1]]

        else: # assume it is an individual index.  Return a singleton dimension.
            fpos, headerlen, _, nsamp, _ = self.index_[idx]
            i0 = fpos + headerlen + self.trace_byteoffset_ * nsamp
            # np.copy is to address unfoc Issue #6: BufferError: cannot close exported pointers exist
            # which occurs with numpy version 1.22.4 and later
            data = np.copy(np.frombuffer(self.mmbxds_, dtype=">i2", offset=i0, count=nsamp))

        return data



    def ct(self, idx):
        """ Return a one or more ct values.  Supports slicing.
        # Get the fourth ct seq/tim value
        one_ct = self.ct(3)
        # Get the fourth and fifth seq/tim values as a list
        two_cts = self.ct(3:5)
        """
        if self.cts_ is None: # Lazily load CTs

            all_cts = tuple(gen_ct(self.fd_.name))
            self.cts_ = [all_cts[idx[4]] for idx in self.index_]
        return self.cts_[idx]



class RadBxdsEx:
    """ Combine one or more RADnh3/RADnh5 channels as specified in the channel spec.
    (currently only supports combining two channels with summation)

    """
    def __init__(self, filename=None, channels=None, stream=None, dtype=None, bxds_class=RadBxds):
        self.rbxds0_ = None
        self.rbxds1_ = None
        self.dtype = None
        self.bxds_class_ = bxds_class
        if filename is not None:
            self.open(filename, channels, stream, dtype)

    def open(self, filename, channels, stream=None, dtype=None):
        """
        filename: bxds file to load
        channels: list of channels to load (a PIK1ChannelSpec object)
        stream: hint at stream type
        """

        # PIK1ChannelSpec(chanout=1, chan0in=1, scalef0=1, chan1in=3, scalef1=1)
        assert channels.scalef0 == 0 or channels.scalef0 == 1
        assert channels.scalef1 == 0 or channels.scalef1 == 1

        if channels.scalef0 == 1 and channels.scalef1 == 1:
            # sum channels
            self.rbxds0_ = self.bxds_class_(filename, channels.chan0in, stream)
            self.rbxds1_ = self.bxds_class_(filename, channels.chan1in, stream)
            self.len_ = min(len(self.rbxds0_), len(self.rbxds1_))
        else:
            # If we're only doing one channel, it better be chan0
            assert channels.scalef0 == 1
            self.rbxds0_ = RadBxds(filename, channels.chan0in, stream)
            self.rbxds1_ = None
            self.len_ = len(self.rbxds0_)

        if dtype is None:
            # Use '>i2' if not combining, otherwise '>i4'
            self.dtype = '>i2' if self.rbxds1_ is None else '>i4'
        else: # Use the user-specified dtype
            self.dtype = dtype

    def __del__(self):
        self.close()

    def close(self):
        self.rbxds0_.close()
        if self.rbxds1_ is not None:
            self.rbxds1_.close()

    def __len__(self):
        # We track len in an internal variable so that we return the right
        # number of CT
        return self.len_

    def __getitem__(self, idx):
        """ Return data for selected radar records combined as an ndarray.

        Records to return can be selected using slice notation.
        See RadBxds.__getitem__ for examples.
        """

        arr1 = self.rbxds0_[idx].astype(self.dtype, copy=False)

        if self.rbxds1_ is None:
            return arr1

        arr2 = self.rbxds1_[idx]

        if arr1.shape != arr2.shape:
            # If they aren't the same size, use the smaller of the two

            # Assume they are the same dimensionality (guaranteed by RadBxds.__getitem__)
            # assert len(arr1.shape) == len(arr2.shape) == 2
            # assume the fast time dimension is the same
            # assert len(arr1.shape) == 1 or arr1.shape[1] == arr2.shape[1]
            if arr1.shape[0] < arr2.shape[0]:
                arr2 = arr2[0:arr1.shape[0], :]
            else:
                arr1 = arr1[0:arr2.shape[0], :]

        arr1 += arr2.astype(self.dtype, copy=False)
        return arr1


    def ct(self, idx):
        """ Return the CT of the first channel """
        if isinstance(idx, slice):
            # Return the shorter of the two sequences
            idx = slice(*idx.indices(self.len_))
        return self.rbxds0_.ct(idx)

class RADjh1Bxds:
    """ Reader for RADjh1 bxds.  This is really just a thin wrapper
    around the numpy.memmap interface for parallelism with RadBxds.
    You might be better off just using np.memmap.
    """

    def __init__(self, filename=None, channel=None, stream=None):
        """ Initialize the reader to access one channel's records """
        self.close()
        if filename is not None:
            self.open(filename, channel, stream)

    def __del__(self):
        self.close()

    def close(self):
        self.cts_ = None
        self.mmbxds_ = None
        self.filename = None
        self.channel_ = None

    def open(self, filename, channel, stream=None):
        """ Open the bxds and load record index for bxds
        filename: bxds filename for low gain channel
        channel: one-based channel number
        stream: (optional) hint for stream type

        """

        assert channel == 1 or channel == 2
        assert os.path.basename(filename) == ('bxds%d' % channel)

        self.channel_ = channel
        self.filename_ = filename

        filesize = os.path.getsize(filename)
        ntraces = filesize // (3200*2) # 3200 samples * 2 bytes per sample

        self.mmbxds_ = np.memmap(filename, dtype='>i2', mode='r', shape=(ntraces, 3200))

    def size(self):
        return self.__len__()

    def __len__(self):
        """ Return the number of traces for this channel """
        return self.mmbxds_.shape[0]

    def __getitem__(self, idx):
        """ View to the underlying memory map  """
        return self.mmbxds_[idx]

    def ct(self, idx):
        """ Return a one or more ct values.  Supports slicing.
        # Get the fourth ct seq/tim value
        one_ct = self.ct(3)
        # Get the fourth and fifth seq/tim values as a list
        two_cts = self.ct(3:5)
        """
        if self.cts_ is None: # Lazily load CTs
            self.cts_ = tuple(gen_ct(self.filename_))

        return self.cts_[idx]


# This could be done but we don't need it.
#class RADjh1BxdsEx(RadBxdsEx):
#    """ Combine one or more RADjh1 channels as specified in the channel spec.
#    (currently only supports combining two channels with summation)
#
#    """
#    def __init__(self, filename=None, channels=None, stream=None, dtype=None, bxds_class=RADjh1Bxds)
#        super(RADjh1BxdsEx, self).__init__(filename, channels, stream, dtype, bxds_class)


def read_1m_gen(bxdsfile, channel, samples_per_trace=3200):
    """
    Set up the generator for reading from a S2_FIL bxdsN.i file
    Expects bxdsN.i to be an array of 2-byte little endian integers,
    typically 3200 samples per trace.
    If the file is not of the expected size, it raises an AssertionError
    Channel is currently required, but hypothetically you could figure it out from the bxdsfile filename.
    """
    nbytes = os.path.getsize(bxdsfile)
    assert nbytes % (2*samples_per_trace) == 0
    ntraces = nbytes // (2*samples_per_trace)

    radargram_in = np.memmap(bxdsfile, dtype='<i2', mode='r', shape=(ntraces, samples_per_trace))

    for ii, arr_trace in enumerate(radargram_in):
        yield Trace(channel=channel, data=arr_trace, ct=CT_t(tim=0, seq=ii))


