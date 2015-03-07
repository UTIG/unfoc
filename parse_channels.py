#!/usr/bin/env python

import collections
# Fields are: output channel, input 1 channel, input 1 scale factor, 
#             input 2 channel, input 2 scale factor.
# function for parsing channel specification
# definition: [input channel
# chanstr = [1,2,3,4,5;6,7,8.0,9,10]
# chanstr = [1,2,3,4,5]
# chanstr = 

PYK1ChannelSpec = collections.namedtuple("PYK1ChannelSpec", ['chanout','chan0in','scalef0','chan1in','scalef1']);

def parse_channels( chanstr ):
    try:
        # Assume it's a simple integer
        chan = int(chanstr)
        return [ PYK1ChannelSpec._make([chan, chan, 1, 0, 0]) ];
        return list_config
    except ValueError as e:
        # If not, that's ok.  
        pass

    # Now parse it as a matlab string
    # Remove spaces and \
    chanstr=chanstr.replace(" ","").replace("\\","");
    # strip beginning/ending brackets
    chanstr=chanstr.lstrip(" [").rstrip(" ];")
    print "chanstr: %r" % (chanstr)
    list_config=[]
    for line in chanstr.split(";"):
        cols=line.split(",")
        if len(cols) != 5:
            raise ValueError("Configuration line ({0:s}) must have 5 values.".format(line))

        # Parse strings into a config.  
        cfgdata = [int(cols[0]), int(cols[1]), float(cols[2]), int(cols[3]), float(cols[4])]
        list_config.append(PYK1ChannelSpec._make(cfgdata))
    return list_config


# Test
def main():
    print parse_channels("1")

    # This should fail
    try:
        parse_channels("1.2")
        exit(1)
    except ValueError as e:
        pass

    print parse_channels("[1,2,3,4,5;6,7,8.0,9,10]")
    print parse_channels("[1,2,3,4,5]")

    for s in "[1,1,1,3,1] [2,2,1,4,1] [5,1,1,0,0;7,3,1,0,0] [6,2,1,0,0;8,4,1,0,0]".split(" "):
        print "s: %s" % s
        print parse_channels(s)

    return 0


if __name__ == "__main__":
    main()

