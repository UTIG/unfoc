#!/usr/bin/env python3

"""
Module for parsing channel specifications.

This is a bit of a relic from the older version of pik1, and a lot of the
functionality doesn't make sense any more.

Its purpose is to map from input data product types to named output channels
to an internal object that the unfocused processor uses to control processing.

Roadmap for improvement plans:
- rename to channels.py
- Remove support for parsing integer channels
- Remove support for linear combinations; just sums?
- Simplify parameters (since scale factor isn't actually supported in current version)
- Ensure that multipol is supported

Basically, the members and meaning of PIK1ChannelSpec are internal and subject to change.

Of course it ain't broken so why fix it?

"""


import logging
import collections
from typing import List, Union



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
                                          'chan1in','scalef1', 'burstnoise_chan0', 'burstnoise_chan1'],
                  defaults=(None,None) );
"""
PIK1ChannelSpec Defines how an output channel is constructed from one or two input channels.

Fields
------
chanout : int
    Output channel number (1-based).
chan0in : int
    First input channel number.
scalef0 : float
    Scale factor for first input channel (usually 1.0).
chan1in : int
    Second input channel number (0 if unused).
scalef1 : float
    Scale factor for second input channel (0 if unused).
burstnoise_chan0 : dict or None
    Burst noise suppression parameters for chan0in (optional).
burstnoise_chan1 : dict or None
    Burst noise suppression parameters for chan1in (optional).
"""

# Default parameters for burst noise on channel 6
burstnoise6 = {
    'median_size': (1, 51), # currently only supports 1d median
    'burst_widths': [15.],
    'detect_thresholds': [25.],
}
# A quick reference for how to set up pik1 to make the UTIG data correctly
# See document multipol_data_format.pdf for more details

UTIG_CHANNELS = {
'MARFA': {

    'LoResInco1': PIK1ChannelSpec(chanout=1, chan0in=1, scalef0=1, chan1in=3, scalef1=1), # sum left and right low gain (XT)
    'LoResInco2': PIK1ChannelSpec(chanout=2, chan0in=2, scalef0=1, chan1in=4, scalef1=1), # sum left and right high gain (XT)

    # 2021-12-21 LoResInco3 and LoResInco4 are being reassigned from difference to summed AT-TX
    #'LoResInco3': PIK1ChannelSpec(chanout=3, chan0in=1, scalef0=1, chan1in=3, scalef1=-1), # diff left and right low gain
    #'LoResInco4': PIK1ChannelSpec(chanout=4, chan0in=2, scalef0=1, chan1in=4, scalef1=-1), # diff left and right high gain
    'LoResInco3': PIK1ChannelSpec(chanout=3, chan0in=1+0x40, scalef0=1, chan1in=3+0x40, scalef1=1), # sum left and right low gain (AT)
    'LoResInco4': PIK1ChannelSpec(chanout=4, chan0in=2+0x40, scalef0=1, chan1in=4+0x40, scalef1=1), # sum left and right high gain (AT)

    # LoResInco5 through LoResInco12 assume an cross-track transmit polarization
    'LoResInco5': PIK1ChannelSpec(chanout=5, chan0in=1, scalef0=1, chan1in=0, scalef1=0), # left low (XT)
    'LoResInco6': PIK1ChannelSpec(chanout=6, chan0in=2, scalef0=1, chan1in=0, scalef1=0), # left high (XT)
    'LoResInco7': PIK1ChannelSpec(chanout=7, chan0in=3, scalef0=1, chan1in=0, scalef1=0), # right low (XT)
    'LoResInco8': PIK1ChannelSpec(chanout=8, chan0in=4, scalef0=1, chan1in=0, scalef1=0), # right high (XT)

    # Added 2021-12-21 for multipol data format
    'LoResInco9': PIK1ChannelSpec(chanout=9, chan0in=5, scalef0=1, chan1in=0, scalef1=0), # cen-fwd low
    'LoResInco10': PIK1ChannelSpec(chanout=10, chan0in=6, scalef0=1, chan1in=0, scalef1=0), # cen-fwd high
    #'LoResInco11': PIK1ChannelSpec(chanout=11, chan0in=7, scalef0=1, chan1in=0, scalef1=0), # cen-aft low
    #'LoResInco12': PIK1ChannelSpec(chanout=12, chan0in=8, scalef0=1, chan1in=0, scalef1=0), # cen-aft high

    # LoResInco13 through LoResInco16 assume an along-track transmit polarization
    'LoResInco13': PIK1ChannelSpec(chanout=13, chan0in=1+0x40, scalef0=1, chan1in=0, scalef1=0), # left low (AT)
    'LoResInco14': PIK1ChannelSpec(chanout=14, chan0in=2+0x40, scalef0=1, chan1in=0, scalef1=0), # left high (AT)
    'LoResInco15': PIK1ChannelSpec(chanout=15, chan0in=3+0x40, scalef0=1, chan1in=0, scalef1=0), # right low (AT)
    'LoResInco16': PIK1ChannelSpec(chanout=16, chan0in=4+0x40, scalef0=1, chan1in=0, scalef1=0), # right high (AT)
    'LoResInco17': PIK1ChannelSpec(chanout=17, chan0in=5+0x40, scalef0=1, chan1in=0, scalef1=0), # cen-fwd low (AT)
    'LoResInco18': PIK1ChannelSpec(chanout=18, chan0in=6+0x40, scalef0=1, chan1in=0, scalef1=0), # cen-fwd high (AT)
    #'LoResInco19': PIK1ChannelSpec(chanout=19, chan0in=7+0x40, scalef0=1, chan1in=0, scalef1=0), # cen-aft low (AT)
    #'LoResInco20': PIK1ChannelSpec(chanout=20, chan0in=8+0x40, scalef0=1, chan1in=0, scalef1=0), # cen-aft high (AT)},
},
'HiCARS2': {
    # For HiCARS, the first two digitizers represent combined data, so output to 1 and 2 not 5 and 6.
    'LoResInco1': PIK1ChannelSpec(chanout=1, chan0in=1, scalef0=1, chan1in=0, scalef1=0), # pass low gain
    'LoResInco2': PIK1ChannelSpec(chanout=2, chan0in=2, scalef0=1, chan1in=0, scalef1=0), # pass high gain
},
}
UTIG_CHANNELS['HERA'] = UTIG_CHANNELS['MARFA']
UTIG_CHANNELS['MPOL'] = UTIG_CHANNELS['MARFA']

def get_utig_channels(chanstr:str, radar:str='MARFA', input_channels:List[int]=None)->PIK1ChannelSpec:
    """ 
    Return a list of PIK1ChannelSpec objects for the given radar and channel names.
    
    Expects a comma-separated list of channels to produce. Case sensitive.
    Radar is a string describing the radar data format.
    Recommend to use unfoc.py autodetect function (get_radar_type)

    If input_channels is provided, this is a list of input channels available in the input
    raw data file. Output channels are filtered according to the
    available channels.
    Note that input_channels and Pik1ChannelSpec both use 1-based indexing.

    Parameters
    ----------
    chanstr : str
        Comma-separated list of channel names, eg "LoResInco1,LoResInco2".
    radar : str, optional
        Radar system name. One of 'MARFA', 'HiCARS2'.
        Default is 'MARFA'.
    input_channels : list of int, optional
        List of available input channels (one-based). If provided,
        channels not present in this list will be excluded.

    Returns
    -------
    list of PIK1ChannelSpec
        Channel configurations matching the requested output channels.

    Raises
    ------
    AssertionError
        If radar is not recognized in UTIG_CHANNELS.

    """
    list_config = [] # type: List[PIK1ChannelSpec]
    assert radar in UTIG_CHANNELS
    for name in chanstr.split(','):
        if name.startswith('LoResInco') and name not in UTIG_CHANNELS[radar]:
            # Silently allow you to specify invalid LoResInco channels
            logging.debug("unfoc::parse_channels: Channel %s isn't known for radar %s", name, radar)
            continue

        p1cs = UTIG_CHANNELS[radar][name]
        if input_channels is not None and (\
           (p1cs.chan0in > 0 and p1cs.chan0in not in input_channels) or \
           (p1cs.chan1in > 0 and p1cs.chan1in not in input_channels)):
                logging.debug("unfoc::parse_channels: skip channel %s", name)
                continue
        # Only add to config if input is available.
        list_config.append(p1cs)
    return list_config

def parse_channels(chanstr):
    """
    Parse a legacy-style channel specification into a list of PIK1ChannelSpec objects.

    Supports two formats:
    - An integer channel number (deprecated).
    - A MATLAB-style semicolon-separated list of 5-column values,
      e.g., "[1,1,1,0,0; 2,2,1,0,0]".

    Parameters
    ----------
    chanstr : int or str
        Integer channel number or formatted string defining multiple channels.

    Returns
    -------
    list of PIK1ChannelSpec
        Parsed list of channel specifications.

    Raises
    ------
    ValueError
        If the string format is invalid (wrong number of values).
    """

    # type: (Union[int, str]) -> List[PIK1ChannelSpec]
    try:
        # Assume it's a simple integer
        chan = int(chanstr)
        logging.warning("Specifying a unfoc channel as a simple integer is deprecated")
        return [PIK1ChannelSpec(chanout=chan, chan0in=chan, scalef0=1, chan1in=0, scalef1=0)];
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
        list_config.append(PIK1ChannelSpec(*cfgdata))
    return list_config
