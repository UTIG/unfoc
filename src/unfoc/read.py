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

from collections import namedtuple, defaultdict, Counter
import itertools
import logging
import os
from pathlib import Path
import struct
import gzip
import mmap
import pickle

import numpy as np

# Pik1ChannelSpec
#import unfoc.parse_channels
from .trace import Trace
from .burst_noise import BurstDenoiser

__version__ = '2.1.0'

radnh3_header = namedtuple('radnh3_header', 'nsamp nchan vr0 vr1 choff ver resvd2')
radnh5_header = namedtuple('radnh5_header', 'nsamp nchan vr0 vr1 choff ver resvd2 absix relix xinc rseq scount tscount')

HEADER_FORMATS = {
    'RADnh3': {
        'header_t': radnh3_header,
        'fmtstr':  '>HBBBBBB',
    },
    'RADnh5': {
        'header_t': radnh5_header,
        'fmtstr': '>HBBBBBBddfLHL',
    }
}






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
CT_full_t = namedtuple('CT_full', 'prj set trn seq year mon day hour min sec hun tim')
def gen_ct(bxdsfile, full=False):
    """ Generate ct data from the ct file associated with bxdsfile
    By default, we only return relevant information. We actually only want the seq and ct tim.
    If full is True, then return all fields.
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

    if full:
        for line in fh:
            fields = line.rstrip().split()
            yield CT_full_t(*fields[0:3], * tuple(map(int, fields[3:])))
    else:
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
    gen = index_RADnhx_bxds(bxdsfile, stream=stream)
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
        - The number of unique channel offsets present in the bxds file
    """
    # mapping of rseq number and the file position of that record
    rseq_to_fpos = defaultdict(list)

    if stream is None:
        stream = get_radar_stream(bxdsfile)

    if stream == 'RADnh3':
        # rseq needs to come from the seq
        genct1 = gen_ct(bxdsfile)
        #f_rseq = lambda ct: ct[0]
        #genct2 = map(f_rseq, genct1)
    else:
        assert stream == 'RADnh5'
        f_rseq = lambda header: header.rseq

    gen1 = index_RADnhx_bxds(bxdsfile, stream=stream, full_header=True)
    gen2 = itertools.islice(gen1, 0, nrecords)
    for ii, (fpos, _, header) in enumerate(gen2):

        if stream == 'RADnh3':
            rseq = next(genct1)[0]
        else:
            # assert stream == 'RADnh5'
            rseq = header.rseq

        rseq_to_fpos[rseq].append((ii, fpos))
        if len(rseq_to_fpos) > 10:
            break

    # Number of channels is the max number of times we see an rseq value
    nchan = max([len(fposlist) for fposlist in rseq_to_fpos.values()])
    # TODO: we can also get this from the nchan field in RADnh3 and RADnh5. assertion check?
    
    for rseq, fposlist in sorted(rseq_to_fpos.items()):
        if len(fposlist) == nchan:
            # Return the first rseq and first position where this was seen.
            return *fposlist[0], rseq, nchan

    raise ValueError("No radar records in %s" % bxdsfile)


def radar_index_summary(bxdsfile:str, stream=None):
    """ Inspect bxds file to determine which records are complete (all digitizers present),
    what type of stream it is (RADnh3 or RADnh5)
    how many digitizers there are
    what the first radar sequence number is
    what the first radar sequence number that has all digitizer records
    what the last radar sequence number is that has all digitizer records
    and whether there are any radar sequence numbers that are missing digitizer records
    and what they are

    Return a dictionary with this information
    """
    p_bxdsfile = Path(bxdsfile)
    stat = p_bxdsfile.stat()
    if stream is None:
        stream = get_radar_stream(bxdsfile)

    rec0, fpos0, rseq0, nchan = sync_radar_start(bxdsfile, stream=stream)

    info = {
        'filename': str(bxdsfile),
        'filesize': stat.st_size,
        'filemtime': int(stat.st_mtime),
        'format': stream,
        'rseq_valid_range': [rseq0, None], # need to figure out the last complete
        'valid_records': 0, # number of records with all radar records
        'nsamples': None,
        'nchan': nchan,
        'incomplete_records': {},
    }
    gen1 = index_RADnhx_bxds(bxdsfile, stream=stream, full_header=True)
    if stream == 'RADnh3':
        ctgen = gen_ct(bxdsfile)
    last_rseq = info['rseq_valid_range'][0]
    #valid_records = 1
    valid_records = 0
    nsamp = None
    rseqs = Counter()
    for ii, (fpos, _, header) in enumerate(gen1):
        if stream == 'RADnh5':
            rseq = header.rseq
        elif stream == 'RADnh3':
            rseq = next(ctgen).seq
        else: # pragma: no cover
            raise ValueError('Processing for %s not yet implemented' % stream)

        rseqs[rseq] += 1
        if rseqs[rseq] == nchan:
            last_rseq = rseq
            valid_records += 1
            del rseqs[rseq]

        if nsamp is None:
            nsamp = header.nsamp
        elif nsamp != header.nsamp:
            logging.warning("Number of samples per trace changed at from %d to %d at rseq=%d",
                            nsamp, header.nsamp, rseq)

    info['rseq_valid_range'][1] = last_rseq
    info['valid_records'] = valid_records
    info['nsamples'] = nsamp

    if rseqs: # if any orphan records (including before or after valid range)
        info['incomplete_records'] = dict(rseqs)

    return info



def index_RADnhx_bxds(input_filename, stream=None, full_header=False, filepos:int=None,
                      buffering:int=-1):
    # type: (str) -> Generator[tuple]
    """ Read the positions of packets within a RADnh3 and RADnh5 bxds file
    and return these as a generator
    Routines can then use this to seek to the correct location in a file.

    if filepos is not None, seek to this file position before indexing

    # TODO: rework this function to always return the full header and make
    all callers grab the member methods
    """

    if stream is None:
        stream = get_radar_stream(input_filename)

    if stream not in HEADER_FORMATS:
        raise ValueError("Invalid stream type %r for file: %s" % (stream, input_filename))

    header_t = HEADER_FORMATS[stream]['header_t']
    fmtstr = HEADER_FORMATS[stream]['fmtstr']
    fmtlen = struct.calcsize(fmtstr)
    # In theory, it's faster to do this here and only compile the format string once.
    header_struct = struct.Struct(fmtstr)

    with open(input_filename, 'rb') as fd:
        if filepos is None:
            fpos = 0
        else:
            fpos = filepos
            fd.seek(fpos)

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
                    assert header.tscount < 100, \
                        "%s fpos=%d header.tscount=%d" % (input_filename, fpos, header.tscount)
                    #tstamps = struct.unpack_from('>{:d}d'.format(tscount), fd.read(8*tscount))
                    fd.seek(8*header.tscount, os.SEEK_CUR)
                    headerlen += 8*header.tscount

            input_samples = header.nsamp
            assert 0 < input_samples < 10000, \
                "%s fpos=%d header.nsamp=%d" % (input_filename, fpos, input_samples)

            # TODO: warn if nsamp != 3200 or 3437?

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



def read_RADnhx_gen(bxds_filename:str, channel:int, stream:str=None, filepos:int=None,
                    buffering:int=-1):
    """ Return a sequence of traces from only one channel, from a bxds
    channel offset is a one-based index.

    Takes a bxds filename as an argument.
    Reads from this file and a ct file alongside it.

    filepos has the same meaning as in index_RADnhx_bxds

    """
    assert 1 <= channel # one-based channel number
    # choff1 is used to calculate the byte offset of radar data within the packet.
    choff1 = (channel - 1) % 2
    # Filter for listed channel offset in radar record
    chfilter = (channel - 1) - choff1

    with open(bxds_filename, 'rb', buffering=buffering) as fd: # for reading traces
        ctgen = gen_ct(bxds_filename)
        radgen = index_RADnhx_bxds(bxds_filename, stream, filepos=filepos, buffering=buffering)
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

class RadBxdsIterator:
    """ Iterator for RadBxds, Radjh1Bxds, RadBxdsEx """
    def __init__(self, bxds):
        self.idx = 0
        self.bxds = bxds

    def __next__(self):
        """ Iterate over the first axis """
        if self.idx >= len(self.bxds):
            raise StopIteration
        value = self.bxds[self.idx]
        self.idx += 1
        return value

class RadBxds:
    """ Reader for random access to traces in a RADnh3 or RADnh5 bxds as numpy arrays
    Future: add support for also auto-detecting RADjh1 bxds
    alternative name ideas:
    
    """
    def __init__(self, filename:str=None, channel:int=None, stream=None, burstnoise=None,
                 validonly=False, indexfile:str=None):
        """ Initialize the reader to access one channel's records
        if burstnoise is a dictionary, pass this to the burst noise calculation
        """
        self.index_ = []
        #self.mmbxds_ = None
        self.fd_ = None
        self.cts_ = None
        self.burstnoise = None
        self.shape = tuple()
        if filename is not None:
            self.open(filename, channel, stream, burstnoise, validonly, indexfile)

    def __del__(self):
        self.close()

    def close(self):
        #if self.mmbxds_ is not None:
        #    self.mmbxds_.close()
        #    self.mmbxds_ = None
        if self.fd_ is not None:
            self.fd_.close()
            self.fd_ = None
        self.cts_ = None
        self.index_ = []


    def open(self, filename:str, channel:int, stream=None, burstnoise=None,
             validonly:bool=False, indexfile:str=None):
        """ Open the bxds and load record index for bxds
        This ignores sequence numbers and assumes that there are no
        out-of-order radar records within a channel.
        In theory we could sort the records by sequence number.

        filename: bxds filename
        channel: one-based radar channel number
        stream: (optional) hint for stream type
        burstnoise: if burstnoise is a dictionary, pass this to the burst noise calculation
        validonly: If validonly is true (will become default in the future), then only return records
        that are valid among all channels

        """

        assert channel >= 1, "Channel parameter should be 1-based"
        self.channel0_ = channel - 1 # zero-based channel index
        self.trace_byteoffset_ = 2 * (self.channel0_ & 1)
        self.fd_ = open(filename, 'rb')
        #self.mmbxds_ = None

        if indexfile is not None:
            self.load_index(indexfile)
        else:
            self.make_index(validonly=validonly, stream=stream)

        # Channel offset that we are expecting to filter for.
        choff = self.channel0_ - (self.channel0_ & 1)
        if len(self.index_[choff]):
            # pick an arbitrary header and assume nsamples is the same for all
            header = self.index_[choff][-1][2]
            self.shape = (len(self.index_[choff]), header.nsamp)
        # otherwise should we set it to empty tuple?


        # Setup burst denoising if in a channel that wants it
        if burstnoise is not None:
            self.burstnoise = BurstDenoiser(**burstnoise)


    def make_index(self, validonly:bool, stream:str=None):
        """ Make index data from data in the file descriptor 
        
        """
        # save all the channels but in a dictionary of lists
        # then we can save this dictionary to the pickle file
        # or perhaps make it a memmap
        # see how this squares against the HDF5 index example
        #@@self.index_ = []
        self.index_ = defaultdict(list)

        #---------------------
        # calculate file size
        prev = self.fd_.tell()
        self.fd_.seek(0, 2)
        filesize = self.fd_.tell()
        self.fd_.seek(prev, 0)
        #---------------------
        filename = self.fd_.name
        if stream is None:
            stream = get_radar_stream(filename)

        if stream == 'RADnh3':
            genct1 = gen_ct(filename)

        if validonly:
            radinfo = radar_index_summary(str(filename), stream=stream)
            #rec0, fpos0, rseq0, nchan = sync_radar_start(str(bxds), stream=s)
            rseq_min, rseq_max = radinfo['rseq_valid_range']
        else:
            rec0, fpos0, rseq0, nchan = 0, 0, None, None
            rseq_min, rseq_max = None, None

        # Channel offset that we are expecting to filter for.
        choff = self.channel0_ - (self.channel0_ & 1)
        # Length of a record's trace data in bytes per sample
        bytes_per_samp = ((self.channel0_ & 1) + 1) * 2
        for ii, item in enumerate(index_RADnhx_bxds(filename, stream=stream, full_header=True)):
            fpos, headerlen, header = item

            if rseq_min is not None: # rseq is specified
                if stream == 'RADnh5':
                    rseq = header.rseq
                elif stream == 'RADnh3':
                    rseq = next(genct1)[0]
                else: #pragma: no cover
                    raise ValueError("Unsupported stream %r" % stream)

                if not (rseq_min <= rseq <= rseq_max):
                    continue


            if header.choff == 0xff: # Replace choff if it is 0xff
                header = header._replace(choff=0x00)

            # Add to index for all channels
            if True: #header.choff == choff:
                lastbyte = fpos + headerlen + bytes_per_samp * header.nsamp
                if lastbyte <= filesize:
                    self.index_[header.choff].append(item + (ii,))

    def save_index(self, indexfile:str):
        """ Save the index data in self.index_ to a file to read later """

        with open(indexfile, 'wb') as fout:
            pickle.dump(dict(self.index_), fout)

    def load_index(self, indexfile:str):
        """ read index data saved with save_index from indexfile into self.index  """

        with open(indexfile, 'rb') as fout:
            self.index_ = pickle.load(fout)
        # invalidate cts() so it has to be reloaded
        self.cts_ = None

    def __getattr__(self, k):
        """ Implement ndarray attributes """
        if k == 'size':
            return self.shape[0] * self.shape[1]
        elif k == 'ndim':
            return len(self.shape)
        elif k == 'dtype':
            return '>i2'
        elif k == 'nbytes':
            return self.shape[0] * self.shape[1] * 2
        elif k == 'T':
            raise NotImplementedError('Transpose view is not yet implemented')
        else:
            raise AttributeError(k)

    def __len__(self):
        """ Return the number of traces for this channel """
        return 0 if len(self.shape) == 0 else self.shape[0]


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
        choff = self.channel0_ - (self.channel0_ & 1)

        if isinstance(idx, slice): # 1D slice case
            indices = idx.indices(len(self.index_[choff]))
            ntraces = len(range(*indices))
            data = np.empty((ntraces, self.index_[choff][0][2].nsamp), dtype=">i2")

            for ii, idxinfo in enumerate(self.index_[choff][idx]):
                fpos, headerlen, header, _ = idxinfo
                i0 = fpos + headerlen + self.trace_byteoffset_ * header.nsamp
                self.fd_.seek(i0, 0)
                data[ii, :] = np.frombuffer(self.fd_.read(header.nsamp<<1), dtype=">i2")
        elif isinstance(idx, tuple):
            # multiple indices -- pass others down to ndarray
            data = self.__getitem__(idx[0])
            if isinstance(idx[0], slice):
                return data[:, idx[1]]
            else: # single integer
                return data[idx[1]]

        else: # assume it is an individual index.  Return a singleton dimension.
            fpos, headerlen, header, _ = self.index_[choff][idx]
            i0 = fpos + headerlen + self.trace_byteoffset_ * header.nsamp
            self.fd_.seek(i0, 0)
            data = np.frombuffer(self.fd_.read(header.nsamp<<1), dtype=">i2")

        if self.burstnoise:
            data = self.burstnoise.denoise(data)

        return data

    def __iter__(self):
        """ Return new iterator """
        return RadBxdsIterator(self)


    def ct(self, idx):
        """ Return a one or more ct values.  Supports slicing.
        # Get the fourth ct seq/tim value
        one_ct = self.ct(3)
        # Get the fourth and fifth seq/tim values as a list
        two_cts = self.ct(3:5)
        """
        if self.cts_ is None: # Lazily load CTs
            choff = self.channel0_ - (self.channel0_ & 1)
            all_cts = tuple(gen_ct(self.fd_.name))
            self.cts_ = [all_cts[idx[3]] for idx in self.index_[choff]]
        return self.cts_[idx]



class RadBxdsEx:
    """ Combine one or more RADnh3/RADnh5 channels as specified in the channel spec.
    (currently only supports combining two channels with summation)

    """
    def __init__(self, filename:str=None, channels=None, stream=None,
                 dtype=None, validonly:bool=True, bxds_class=RadBxds,
                 indexfile:str=None):
        self.rbxds0_ = None
        self.rbxds1_ = None
        self.dtype = None
        self.bxds_class_ = bxds_class

        if filename is not None:
            self.open(filename, channels, stream, dtype, validonly, indexfile)

    def open(self, filename:str, channels, stream=None, dtype=None, validonly:bool=True,
             indexfile:str=None):
        """
        filename: bxds file to load
        channels: list of channels to load (a PIK1ChannelSpec object)
        stream: hint at stream type
        indexfile: if provided, load the index from this file
        """

        # PIK1ChannelSpec(chanout=1, chan0in=1, scalef0=1, chan1in=3, scalef1=1)
        assert channels.scalef0 == 0 or channels.scalef0 == 1
        assert channels.scalef1 == 0 or channels.scalef1 == 1

        if channels.scalef0 == 1 and channels.scalef1 == 1:
            # sum channels
            self.rbxds0_ = self.bxds_class_(filename, channels.chan0in, stream, burstnoise=channels.burstnoise_chan0, validonly=validonly, indexfile=indexfile)
            self.rbxds1_ = self.bxds_class_(filename, channels.chan1in, stream, burstnoise=channels.burstnoise_chan1, validonly=validonly, indexfile=indexfile)
            assert self.rbxds0_.shape == self.rbxds1_.shape, "Radar data for %s doesn't have same shape." % (filename)
        else:
            # If we're only doing one channel, it better be chan0
            assert channels.scalef0 == 1
            self.rbxds0_ = self.bxds_class_(filename, channels.chan0in, stream, burstnoise=channels.burstnoise_chan0, validonly=validonly, indexfile=indexfile)
            self.rbxds1_ = None

        if dtype is None:
            # Use '>i2' if not combining, otherwise '>i4'
            self.dtype = '>i2' if self.rbxds1_ is None else '>i4'
        else: # Use the user-specified dtype
            self.dtype = dtype

    def save_index(self, indexfile:str):
        """ Save index for each file using underlying index file.
        Since open() can only be called with one input filename,
        indexfile is also only one filename
        """
        self.rbxds0_.save_index(indexfile)

        if self.rbxds1_ is not None:
            self.rbxds1_.save_index(indexfile)

    def load_index(self, indexfile:str):
        self.rbxds0_.load_index(indexfile)

        if self.rbxds1_ is not None:
            self.rbxds1_.load_index(indexfile)

    def __del__(self):
        self.close()

    def close(self):
        if self.rbxds0_ is not None:
            self.rbxds0_.close()
        if self.rbxds1_ is not None:
            self.rbxds1_.close()

    def __len__(self):
        # We track len in an internal variable so that we return the right
        # number of CT
        return self.rbxds0_.shape[0]

    def __getattr__(self, k):
        """ Implement ndarray attributes """
        if k in ('size', 'shape', 'ndim'):
            return getattr(self.rbxds0_, k)
        elif k == 'nbytes':
            if isinstance(self.dtype, str):
                element_size = np.dtype(self.dtype).itemsize
            else:
                element_size = self.dtype.itemsize
            return self.rbxds0_.shape[0] * self.rbxds0_.shape[1] * element_size
        elif k == 'T':
            raise NotImplementedError('Transpose view is not yet implemented')
        else:
            raise AttributeError(k)


    def __getitem__(self, idx):
        """ Return data for selected radar records combined as an ndarray.

        Records to return can be selected using slice notation.
        See RadBxds.__getitem__ for examples.
        """

        arr1 = self.rbxds0_[idx].astype(self.dtype, copy=False)

        if self.rbxds1_ is None:
            return arr1

        arr2 = self.rbxds1_[idx]
        assert arr1.shape == arr2.shape
        '''
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
        '''
        arr1 += arr2.astype(self.dtype, copy=False)
        return arr1

    def __iter__(self):
        """ Return new iterator """
        return RadBxdsIterator(self)



    def ct(self, idx):
        """ Return the CT of the first channel """
        if isinstance(idx, slice):
            # Return the shorter of the two sequences
            idx = slice(*idx.indices(self.rbxds0_.shape[0]))
        return self.rbxds0_.ct(idx)

class RADjh1Bxds:
    """ Reader for RADjh1 bxds.  This is really just a thin wrapper
    around the numpy.memmap interface for parallelism with RadBxds.
    You might be better off just using np.memmap.
    """

    def __init__(self, filename=None, channel=None, stream=None, indexfile:str=None):
        """ Initialize the reader to access one channel's records
        indexfile is not used but included for interface compatibility
        with other *Bxds classes
        """
        self.fd_ = None
        self.index_ = [] # for compatibility with RadBxds
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
        if self.fd_ is not None:
            self.fd_.close()
            self.fd_ = None

    def open(self, filename:str, channel:int, stream=None):
        """ Open the bxds and load record index for bxds
        filename: bxds filename for low gain channel
        channel: one-based channel number
        stream: (optional) hint for stream type

        """

        assert channel in (1, 2), "RADjh1 only has channels 1 and 2"
        assert os.path.basename(filename) == ('bxds%d' % channel), \
               "Mismatch between filename and channel"

        self.channel_ = channel

        self.filename_ = filename

        filesize = os.path.getsize(filename)
        assert filesize % (3200*2) == 0, "File isn't expected size: %d bytes" % (filesize)
        ntraces = filesize // (3200*2) # 3200 samples * 2 bytes per sample

        #self.fd_ = open(filename, 'rb')
        self.mmbxds_ = np.memmap(filename, dtype='<i2', mode='r', shape=(ntraces, 3200))

    def __getattr__(self, k:str):
        """ Implement ndarray attributes """
        if k == 'size':
            return self.mmbxds_.shape[0] * self.mmbxds_.shape[1]
        elif k == 'ndim':
            return len(self.mmbxds_.shape)
        elif k == 'dtype':
            return '<i2'
        elif k == 'nbytes':
            return self.mmbxds_.shape[0] * self.mmbxds_.shape[1] * 2
        elif k == 'shape':
            return self.mmbxds_.shape
        elif k == 'T':
            raise NotImplementedError('Transpose view is not yet implemented')
        else:
            raise AttributeError(k)

    def __iter__(self):
        """ Return new iterator """
        return RadBxdsIterator(self)


    def __len__(self):
        """ Return the number of traces for this channel """
        return self.mmbxds_.shape[0]

    def __getitem__(self, idx):
        """ View to the underlying memory map  """
        return self.mmbxds_[idx]

    def save_index(self, indexfile:str):
        """ placeholder so tests pass.  There is
        no index to save or load. touch the file? """
        pass
        #with open(indexfile, 'wt') as fout:
        #    pass

    def load_index(self, indexfile:str):
        """ placeholder so tests pass.  There is nothing to load """
        pass
        #with open(indexfile, 'rt') as fin:
        #    pass
        # Don't forget to invalidate self.cts_ if anything gets added here

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
    # TODO: do iteration

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


