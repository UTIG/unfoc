#!/usr/bin/env python3


"""
radbxds.py

Read a radar bxds file and return traces
Read the structure of a radar file
"""

from collections import namedtuple
import itertools
import logging
import os
import struct
import gzip
import mmap

import numpy as np




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




def get_radar_stream(filename):
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




def get_radar_type(bxdsfile, nrecords=1000, stream=None):
    """ Inspect a raw datafile and detect what type of radar it came from
    returns 'HiCARS2' if it is a 1-antenna radar, or 'MARFA' if it is a 2-antenna
    radar.
    """
    choffs = {}
    gen = index_RADnhx_bxds_mmap(bxdsfile, stream=stream)
    # fpos, headerlen, header.choff, header.nsamp
    for ii, (_, _, choff, _) in enumerate(gen):
        choffs[choff] = choffs.get(choff, 0) + 1
        if ii >= nrecords or len(choffs) >= 3:
            break
    gen.close()

    if len(choffs) >= 2:
        return 'MARFA'
    else:
        assert len(choffs) == 1
        return 'HiCARS2'




# Read individual traces out of RADnh3 or RADnh5 file
# TODO: This doesn't yet filter on channels, which should be OK - the
# downstream users ALSO check channel.

def index_RADnhx_bxds(input_filename, stream=None):
    # type: (str) -> Generator[tuple]
    """ Read the positions of packets within a RADnh3 and RADnh5 bxds file
    and return these as a generator
    Routines can then use this to seek to the correct location in a file.
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

            yield fpos, headerlen, header.choff, header.nsamp #, ctdata)
            # Skip the rest of this record without yielding any values
            fd.seek(4*input_samples, os.SEEK_CUR)



# Read individual traces out of RADnh3 or RADnh5 file
# TODO: This doesn't yet filter on channels, which should be OK - the
# downstream users ALSO check channel.
def index_RADnhx_bxds_mmap(input_filename, stream=None):
    # type: (str) -> Generator[tuple]
    """ Read the positions of packets within a RADnh3 and RADnh5 bxds file
    and return these as a generator
    Routines can then use this to seek to the correct location in a file.
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
    choff1 = (channel - 1) % 2
    # Filter for listed channel offset in radar record
    chfilter = (channel - 1) - choff1

    try:
        fd = open(bxds_filename, 'rb')
        mmbxds = mmap.mmap(fd.fileno(), 0, access=mmap.ACCESS_READ)

        ctgen = gen_ct(bxds_filename)

        for radinfo, ctinfo in zip(index_RADnhx_bxds_mmap(bxds_filename, stream), ctgen):
            fpos, headerlen, rchoff, nsamp = radinfo
            if rchoff == 0xff:
                rchoff = 0

            if chfilter != rchoff:
                continue

            i0 = fpos + headerlen + choff1 * nsamp * 2
            try:
                trace1 = np.frombuffer(mmbxds, dtype='>i2', offset=i0, count=nsamp)
            except ValueError as e:
                break # Not enough data

            yield Trace(channel, trace1, ctinfo)
    finally:
        mmbxds.close()
        fd.close()


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


    def open(self, filename, channel, stream=None, do_ct=False):
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
        choff = self.channel0_ - self.channel0_ & 1
        # Length of a record's trace data in bytes per sample
        bytes_per_samp = ((self.channel0_ & 1) + 1) * 2

        for ii, item in enumerate(index_RADnhx_bxds_mmap(filename, stream=stream)):
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
        """

        if isinstance(idx, slice):
            indices = idx.indices(len(self.index_))
            ntraces = len(range(*indices))
            data = np.empty((ntraces, self.index_[0][3]), dtype=">i2")

            for ii, idxinfo in enumerate(self.index_[idx]):
                fpos, headerlen, _, nsamp, _ = idxinfo
                i0 = fpos + headerlen + self.trace_byteoffset_ * nsamp
                data[ii, :] = np.frombuffer(self.mmbxds_, dtype=">i2", offset=i0, count=nsamp)

        else: # assume it is an individual index.  Return a singleton dimension.
            fpos, headerlen, _, nsamp, _ = self.index_[idx]
            i0 = fpos + headerlen + self.trace_byteoffset_ * nsamp
            data = np.frombuffer(self.mmbxds_, dtype=">i2", offset=i0, count=nsamp)

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

class RADjh1Bxds:
    """ Reader for RADjh1 bxds.
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

