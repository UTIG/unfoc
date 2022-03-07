#!/usr/bin/env python3

import logging
import collections
try:
    import typing
    from typing import List, Union
except ImportError: #pragma: no cover
    pass # Not installed on melt ...



# This is a way of mapping from input channels to output channels.
# An output channel is a linear combination of any two input channels.
# * chanout (int) - output channel number
# * chan0in (int) - first input channel number
# * scalef0 (float) - first input channel weight
# * chan1in (int) - second input channel number
# * scalef1 (float) - second input channel weight


# We use a NamedTuple because it's more lightweight than a class, and
# because we aren't encoding any behaviors (methods)
# A namedtuple is also immutable so you can't accidentally change the data
PIK1ChannelSpec = collections.namedtuple("PIK1ChannelSpec",
                                         ['chanout','chan0in','scalef0',
                                          'chan1in','scalef1']);

# A quick reference for how to set up pik1 to make the UTIG data correctly

UTIG_CHANNELS = {
'MARFA': {
    'LoResInco1': PIK1ChannelSpec(chanout=1, chan0in=1, scalef0=1, chan1in=3, scalef1=1), # sum left and right low gain
    'LoResInco2': PIK1ChannelSpec(chanout=2, chan0in=2, scalef0=1, chan1in=4, scalef1=1), # sum left and right high gain
    'LoResInco3': PIK1ChannelSpec(chanout=3, chan0in=1, scalef0=1, chan1in=3, scalef1=-1), # diff left and right low gain
    'LoResInco4': PIK1ChannelSpec(chanout=4, chan0in=2, scalef0=1, chan1in=4, scalef1=-1), # diff left and right high gain
    'LoResInco5': PIK1ChannelSpec(chanout=5, chan0in=1, scalef0=1, chan1in=0, scalef1=0), # left low
    'LoResInco6': PIK1ChannelSpec(chanout=6, chan0in=2, scalef0=1, chan1in=0, scalef1=0), # left high
    'LoResInco7': PIK1ChannelSpec(chanout=7, chan0in=3, scalef0=1, chan1in=0, scalef1=0), # right low
    'LoResInco8': PIK1ChannelSpec(chanout=8, chan0in=4, scalef0=1, chan1in=0, scalef1=0), # right high
},
'HiCARS2': {
    'LoResInco1': PIK1ChannelSpec(chanout=1, chan0in=1, scalef0=1, chan1in=0, scalef1=0), # pass low gain
    'LoResInco2': PIK1ChannelSpec(chanout=2, chan0in=2, scalef0=1, chan1in=0, scalef1=0), # pass high gain
},
}
UTIG_CHANNELS['HERA'] = UTIG_CHANNELS['MARFA']


def get_utig_channels(chanstr, radar='MARFA'):
    """ Expects a comma-separated list of channels to produce. Case sensitive.
    Radar is a string describing the radar data format.
    Recommend to use unfoc.py autodetect function (get_radar_type)
    """
    # type: (str) -> List[PIK1ChannelSpec]
    list_config = [] # type: List[PIK1ChannelSpec]
    for name in chanstr.split(','):
        list_config.append(UTIG_CHANNELS[radar][name])
    return list_config

def parse_channels(chanstr):
    # type: (Union[int, str]) -> List[PIK1ChannelSpec]
    try:
        # Assume it's a simple integer
        chan = int(chanstr)
        logging.warning("Specifying a unfoc channel as a simple integer is deprecated")
        return [PIK1ChannelSpec._make([chan, chan, 1, 0, 0])];
    except ValueError as e:
        # If not, that's ok.
        pass

    # Now parse it as a matlab string
    # Remove spaces and \
    chanstr = chanstr.replace(" ","").replace("\\","");
    # strip beginning/ending brackets
    chanstr = chanstr.lstrip(" [").rstrip(" ];")
    list_config = [] # type: List[PIK1ChannelSpec]
    for line in chanstr.split(";"):
        cols = line.split(",")
        if len(cols) != 5:
            raise ValueError("Configuration line ({0:s}) must have 5 values.".format(line))

        # Parse strings into a config.
        cfgdata = [int(cols[0]), int(cols[1]), float(cols[2]), int(cols[3]), float(cols[4])]
        list_config.append(PIK1ChannelSpec._make(cfgdata))
    return list_config
