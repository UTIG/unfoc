#!/usr/bin/env python3

""" Read a RADnh3 or RADnh5 bxds file, and export
binary files for each channel's radar traces and metadata.
"""

import os
import sys
from pathlib import Path
import struct
import datetime
import logging

import numpy as np
import h5py

sys.path.insert(1, '../src')
import unfoc

def main():
    WAIS = os.getenv('WAIS')
    #bxdsfile = Path(WAIS) / 'orig/xlob/TAM/MKB2o/BDMG01a/RADnh5/bxds'
    bxdsfile = Path(WAIS) / 'orig/xlob/THW2/UBH0c/Y34b/RADnh5/bxds'
    outfile = Path('out1') / 'radar.h5'

    logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)

    export1(bxdsfile, outfile)


def export1(bxdsfile, outfile):
    """ A simple export scheme that duplicates some data and
    throws away some internal timing metadata """
    #radsnm, channels = unfoc.get_radar_type(bxdsfile)


    # Determine file layout
    idx_summary = unfoc.radar_index_summary(bxdsfile)
    print(idx_summary)

    Path(outfile).parent.mkdir(exist_ok=True, parents=True)
    with h5py.File(outfile, 'w') as f:
        # Allocate radargram
        # Allocate datasets for ct and ELSA metadata

        # Write attributes for PST
        write_attributes(f, idx_summary)

        radar_index = write_metadata_radnh5(f, bxdsfile, idx_summary)

        write_radar(f, bxdsfile, radar_index)




def write_radar(f, bxdsfile:Path, radar_index, chunksize:int=1000, limit:int=None):
    logging.debug("write_radar from %s", str(bxdsfile))
    shape = (f.attrs['number_of_traces'],
             f.attrs['number_of_digitizers'],
             f.attrs['channels_per_digitizer'],
             f.attrs['samples_per_trace'])
    f.create_dataset('traces', shape=shape,
                     dtype='>i2',
                     chunks=(chunksize, shape[1], shape[2], shape[3]),
                     # Compress with gzip, but don't try that hard
                     #compression='gzip', compression_opts=1,
                     )

    f['traces'].description = 'Radar digitizer sample data'

    # Rtime is statistics about an array of timestamps of each constituent record in this stack,
    # interval mean and std deviation
    # Convert to 32-bit float because why would you actually need all that precision?
    # We probably should have written it this way in the first place
    metadata_rtime = np.zeros(shape=shape[0:2] + (4,), dtype='>f')


    write_chunks = f['traces'].iter_chunks()
    #for s in write_chunks:
    nsamples = shape[-1]
    with open(bxdsfile, 'rb') as fhbxds:
        for chunknum, s in enumerate(write_chunks):
            logging.debug("chunk %d %r", chunknum, s)
            chunk = np.empty(f['traces'][s].shape, dtype='>i2')

            for ii, tracenum in enumerate(range(s[0].start, s[0].stop, s[0].step)):
                for jj, digitizernum in enumerate(range(s[1].start, s[1].stop, s[1].step)):
                    fpos, valid = radar_index[tracenum, digitizernum]

                    if not valid:
                        chunk[ii, jj, ...] = 0
                        metadata_rtime[tracenum, digitizernum, :] = 0.
                        continue

                    fhbxds.seek(int(fpos) - 4*8) # back up so we can get the rtimes
                    rtimes = np.fromfile(fhbxds, dtype='>d', count=4)
                    metadata_rtime[tracenum, digitizernum, :] = rtimes
                    traces = np.fromfile(fhbxds, dtype='>i2', count=2*nsamples).reshape((2, nsamples))
                    chunk[ii, jj, :, :] = traces

            f['traces'][s] = chunk

            if limit is not None and chunknum >= limit:
                break

    # Write the rtimes arrays
    f.create_dataset('stack_timing/std_dev', data=metadata_rtime[:, :, 0:2], compression='gzip')
    f['stack_timing/std_dev'].attrs['description'] = 'Standard deviation of interval between timestamps of traces used stack, in seconds'
    f.create_dataset('stack_timing/mean', data=metadata_rtime[:, :, 2:4], compression='gzip')
    f['stack_timing/mean'].attrs['description'] = 'Mean of interval between timestamps of traces used in stack, in seconds'

    # Write the validity array
    traces_valid = np.zeros(shape[0], dtype='u1')
    for digitizernum in range(shape[1]):
        traces_valid += radar_index['valid'][:, digitizernum] << digitizernum
    f.create_dataset('traces_valid', data=traces_valid, compression='gzip')
    f['traces_valid'].attrs['description'] = \
        'Each element is a bit vector indicating digitizer data are valid. ' \
        'The least significant bit corresponds to the first digitizer. ' \
        'A 1 bit in the position indicates that this digitizer record is valid.'




def write_metadata_radnh5(f, bxdsfile:Path, radar_idx_summary):
    """ Write the metadata contained in bxdsfile to the hdf5 file handle """
    logging.debug("write_metadata_radnh5")
    gen1 = unfoc.index_RADnhx_bxds(bxdsfile, full_header=True)
    gen2 = unfoc.gen_ct(str(bxdsfile), full=True)
    # Allocate metadata arrays
    myshape = (f.attrs['number_of_traces'], f.attrs['number_of_digitizers'])

    choff_notes = """
choff[4:0] encodes the channel offset, a number between 0 and 31.  Bits choff[7:5] are reserved
for flags relating to channel recording.  choff[5]=1 indicates that the data was recorded using
more than one set of digitizers.

If choff[6]=0, then the radar system was acquired with the primary transmit polarization.
If choff[6]=1, then the radar system acquired this data with the secondary transmit polarity
stimulus. Downstream data processing may elect to select a subset of records based on this bit.

nchan indicates the number of channels in this record. Typically, nchan=2 if the digitizer used is the PXIe-5122.
If nchan=2 and choff[4:0]=0, then channels 0 and 1 are contained in the record.
If nchan=2 and choff[4:0]=2, then channels 2 and 3 are contained in the record.
If nchan=2 and choff[4:0]=4, then channels 4 and 5 are contained in the record.

By convention, channels are ordered from left to right, then low gain to high gain.

In a two channel system, the order is: [0] low gain, [1] high gain

In a four channel system:
[0] low gain left, [1] high gain left,
[2] low gain right, [3] high gain right

In a six channel system:
[0] low gain left, [1] high gain left,
[2] low gain right, [3] high gain right
[4] low gain center, [5] high gain center

To maintain backward compatibility with previous RADnh3 data, the value 0xff in
the choff field indicates that exactly two channels (generally on one digitizer
card) are being digitized by this acquisition system. I have verified that the
current acquisition code will set this field appropriately any time this is true.
"""

    # TODO: 2-tuple at beginning
    radnh5_dtype = [
        (('Number of samples per trace', 'nsamp'), 'u2'),
        (('Number of channels in this digitizer record', 'nchan'), 'u1'),
        (('Channel offset for this record' + choff_notes, 'choff'), 'u1'),
        (('Voltage range flag for channel 0 (usually low gain) in this digitizer record', 'vr0'), 'u1'),
        (('Voltage range flag for channel 1 (usually high gain) in this digitizer record', 'vr1'), 'u1'),
        (('Packet format - For RADnh3 packets, this byte is 0.  For RADnh5, it is 5.', 'ver'), 'u1'),
        (('Reserved byte', 'resvd2'), 'u1'),
        (('absoluteInitialX is the timestamp in seconds of the first fetched sample. '
          'This timestamp is comparable between records and acquisitions', 'absix'), 'd'),
        (('relativeInitialX is the time in seconds from the trigger to the '
          'first sample in the acquired waveform', 'relix'), 'd'),
        (('xIncrement indicates the time in seconds between two samples in the acquired waveform', 'xinc'), 'f'),
        (('Digitizer radar record sequence number of the first trace used in this stack', 'rseq'), 'u4'),
        (('Stack count - Number of onboard stacks incorporated into this digitizer record', 'scount'), 'u2'),
        (('Timestamp count (s)', 'tscount'), 'u4'),
    ]

    ct_dtype = [
        (('ELSA sequence number', 'elsa__seq'), 'u4'),
        (('ELSA clk year', 'elsa__clk__year'), 'u2'),
        (('ELSA clk month', 'elsa__clk__mon'), 'u1'),
        (('ELSA clk day of month', 'elsa__clk__day'), 'u1'),
        (('ELSA clk hour', 'elsa__clk__hour'), 'u1'),
        (('ELSA clk minute', 'elsa__clk__minute'), 'u1'),
        (('ELSA clk seconds', 'elsa__clk__seconds'), 'u1'),
        (('ELSA clk hundredths of a second', 'elsa__clk__hun'), 'u1'),
        (('ELSA counter-timer timestamp (1/100_000 seconds)', 'elsa__tim'), 'u4'),
    ]

    radar_index_dtype = [
        ('fpos', 'u8'), # file position of the radar samples payload
        ('valid', 'u1'), # Is this entry valid?
    ]

    metadata = np.zeros(shape=myshape, dtype=radnh5_dtype)
    metadata_ct = np.zeros(shape=myshape, dtype=ct_dtype)
    radar_index = np.zeros(shape=myshape, dtype=radar_index_dtype)
    psts = set()

    logging.debug("write_metadata_radnh5: iterating over metadata")
    rseq0, rseq1 = radar_idx_summary['rseq_valid_range']
    for ii, ((fpos, headerlen, header), ctrec) in enumerate(zip(gen1, gen2)):
        if not (rseq0 <= header.rseq <= rseq1):
            continue
        assert header.scount == 32, "This code will break if stack count != 32"
        assert header.tscount == 4, "The rtimes code will break if this is not 4"
        idx0 = (header.rseq - rseq0) // 32
        idx1 = header.choff // 2
        metadata[idx0, idx1] = header

        metadata_ct[idx0, idx1] = ctrec[3:]

        radar_index[idx0, idx1] = (fpos + headerlen, 1)

    for (title, name), _ in ct_dtype:
        name1 = name.replace('__', '/') # put in a subgroup
        logging.debug("writing %s (%s)", name1, title)
        f.create_dataset(name1, data=metadata_ct[name], compression='gzip')
        f[name1].attrs['description'] = title


    for (title, name), _ in radnh5_dtype:
        if name == 'resvd2':
            continue
        logging.debug("writing %s (%s)", name, title)
        f.create_dataset(name, data=metadata[name], compression='gzip')
        f[name].attrs['description'] = title

    return radar_index



def write_attributes(f, idx_summary):
    """ write the attributes to the hdf5 file root dataset """
    f.attrs['stream_name'] = idx_summary['format']
    f.attrs['source_filesize'] = idx_summary['filesize']
    f.attrs['source_modified_time'] = datetime.datetime.fromtimestamp(idx_summary['filemtime']).strftime("%c")
    f.attrs['number_of_traces'] = idx_summary['valid_records']
    f.attrs['number_of_digitizers'] = idx_summary['nchan']
    f.attrs['channels_per_digitizer'] = 2
    # should we also just make a shape attribute or nah?
    f.attrs['samples_per_trace'] = idx_summary['nsamples']
    # incomplete records?

if __name__ == "__main__":
    main()

