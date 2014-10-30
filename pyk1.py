#!/bin/env python2.7
#
# Dechirp the JPL and MIRS radars.
#
# All output files are 4-byte network-order.
#

import sys
import argparse
import struct
import numpy
from scipy import signal
# Disable for production work (fewer dependencies)
#import matplotlib.pyplot as plt

################################################
# Enable metadata index writing.  
# This should normally be disabled only for legacy testing.
# It does not make things slower to output this data.
enable_meta_index=True
################################################

def cinterp(Data, index):
    r = (numpy.abs(Data[index-1]) + numpy.abs(Data[index+1])) / 2
    t1 = numpy.angle(Data[index-1])
    t2 = numpy.angle(Data[index+1])
    if (numpy.abs(t1 - t2) > numpy.pi):
        t1 = t1 + 2 * numpy.pi
    theta = (t1 + t2) / 2
    Data[index] = r * (numpy.cos(theta) + 1j * numpy.sin(theta))
    return Data


def denoise_and_dechirp(Stacked, Rchirp, blanking, truncSweepLength):
    Stacked[0:blanking] = numpy.zeros(blanking)

#find peak energy below blanking samples
## [n,m]=sort(Stacked);
## shifter=abs((m(truncSweepLength)));
    shifter=int(numpy.median(numpy.argmax(Stacked)));
    Stacked=numpy.roll(Stacked,-shifter);
        
    DFT = numpy.fft.fft(signal.detrend(Stacked))
    # Remove five samples per cycle problem
    DFT = cinterp(DFT, truncSweepLength * (1.0/5))
    DFT = cinterp(DFT, truncSweepLength * (1 - 1.0/5))
    # Remove the first harmonic for five samples
    DFT = cinterp(DFT, truncSweepLength * (2.0/5))
    DFT = cinterp(DFT, truncSweepLength * (1 - 2.0/5))

    # Do the dechirp
    Product = numpy.multiply(Rchirp, DFT)
    Dechirped = numpy.fft.ifft(Product)
    Dechirped = numpy.roll(Dechirped,shifter)
    return Dechirped

def read_and_stack_RADnh3(InputName,ChannelNum,SweepLength,StackDepth=10):
    fd = open(InputName, 'r')
    stacked = numpy.zeros(SweepLength)
    while True:
        for i in range(0,StackDepth):
# Skip "header"
            numpy.fromfile(fd, dtype='>i4', count=2)
            if (ChannelNum == 1):
                stack = numpy.fromfile(fd, dtype='>i2', count=SweepLength)
# Skip other channel
            numpy.fromfile(fd, dtype='>i2', count=SweepLength)
            if (ChannelNum == 2):
                stack = numpy.fromfile(fd, dtype='>i2', count=SweepLength)
            if (stack.size != SweepLength):
                return
            stacked = stacked + stack
        yield stacked/StackDepth
        stacked = numpy.zeros(SweepLength)


def main(argv):
    parser = argparse.ArgumentParser(description='Pulse compress radar data')
# filename of 2-byte radar file
    parser.add_argument('--InputName', required=True)
# Channel Number to compress (counts from 1)
    parser.add_argument('--ChannelNum', default=0, type=int)
# Length of each sweep (in samples)
    parser.add_argument('--SweepLength', required=True, type=int)
# Length of each output sweep (in samples)
    parser.add_argument('--truncSweepLength', required=True, type=int)
# filename of real part of raw compressed waveform
    parser.add_argument('--IName')
# filename of imaginary part of raw compressed waveform
    parser.add_argument('--QName')
# filename of recitified output
    parser.add_argument('--MagName')
# filename of phase output
    parser.add_argument('--PhsName')
# filename for meta data output
    parser.add_argument('--MetaName')
# coherent stacking depth for this output
    parser.add_argument('--StackDepth', required=True, type=int)
# incoherent stacking depth for this output
    parser.add_argument('--IncoDepth', required=True, type=int)
# stacks are centered on (or before) multiples of CenterMult
    parser.add_argument('--CenterMult', required=True, type=int)
# Only output a sweep if MaxDepth sweeps could be stacked.
    parser.add_argument('--MaxDepth', required=True, type=int)
# Output only sweeps with centers on or after this sweep
    parser.add_argument('--StartSweep', type=int, default=1)
# Output only sweeps with centers before or on this sweep
    parser.add_argument('--EndSweep', type=int)
# Output sweep samples starting with this one
    parser.add_argument('--StartSamp', type=int, default=1)
# Output sweep samples ending with this one
    parser.add_argument('--EndSamp', type=int)
# Output scale default is 1000*dB
    parser.add_argument('--Scale', type=int, default=20000)
# Samples at the top of the record to blank out
    parser.add_argument('--blanking', type=int, default=50)
# Stream name
    parser.add_argument('--StreamName', required=True)

    args = parser.parse_args()

    if args.EndSamp or args.EndSamp == 0 > args.SweepLength:
        args.EndSamp = args.SweepLength
    if args.StartSamp > args.EndSamp:
        sys.exit('pyk1: bad start and end samples.')

    if args.StackDepth <= 0 or args.CenterMult <= 0 or args.MaxDepth <= 0:
        sys.exit('pyk1: all depths must be positive')

    if args.SweepLength <= 0:
        sys.exit('pyk1: sweep length must be positive')

    record_increment=args.StackDepth*args.IncoDepth
    record=record_increment/2

# Obtain reference chirp
    I = numpy.array([-63, -92, -109, -75, -87, -50, -116, -154, -22, 68, 141, -610, 1461, 3807, -6147, -5375, 10651, -4412, -9810, 15386, -3070, -14499, 15130, 3677, -15935, 3743, 13362, -5884, -13301, 8455, 12542, -8744, -11977, 5105, 13754, -961, -14342, -5184, 10294, 12194, -2709, -14352, -8807, 5965, 15350, 8368, -6605, -14990, -11515, 196, 11276, 15490, 10300, -645, -10730, -15307, -13379, -7342, 377, 7264, 11662, 13435, 11530, 3243, -4865, -4427, -3233, -4000, -2472, -2498, -2361, -1230, -1311, -618, -578, -569, -121, -319, 206, 328, 436, 613, 318, 514, 353, 277, 221, 34, 250, 132, 199, 189, 75, 190, 65, 106, 19, -64, -14, -117])
    rchirp = numpy.flipud(I)
    Rchirp = numpy.fft.fft(rchirp,n=args.truncSweepLength)

# Compute Hamming Filter
    From = round(2.5 * args.truncSweepLength / 50)
    To = round(17.5 * args.truncSweepLength / 50)
    hamming = numpy.sin(numpy.linspace(0,1,num=To-From+1) * numpy.pi)
    filter = numpy.hstack((numpy.zeros(From), hamming*2, numpy.zeros(args.truncSweepLength-2*To-1), numpy.zeros(hamming.size), numpy.zeros(From-1)))

## Filter Chirp
    # disable for production work
    # plt.plot(filter, range(0,filter.size))
    Rchirp = numpy.multiply(Rchirp, filter)

## Open Input and Output files

    MetaOutFD = open(args.MetaName, 'w')
    MetaOutFD.write('#InputName = "' + args.InputName + '"\n')

    MagOutFD = open(args.MagName, 'w')
    MetaOutFD.write('#MagName = "' + args.MagName + '"\n')

    if args.PhsName:
        PhsOutFD = open(args.PhsName, 'w')
        MetaOutFD.write('#PhsName = "' + args.PhsName + '"\n')

    ###### FIXME? - request hdr file from xlob and read/forward
    ###### Make this Meta similar
    MetaOutFD.write("#ChannelNum = " + str(args.ChannelNum) + "\n")
    MetaOutFD.write("#StackDepth = " + str(args.StackDepth) + "\n")
    MetaOutFD.write("#IncoDepth = " + str(args.IncoDepth) + "\n")
    MetaOutFD.write("#CenterMult = " + str(args.CenterMult) + "\n")
    MetaOutFD.write("#MaxDepth = " + str(args.MaxDepth) + "\n")
    MetaOutFD.write("#StartSweep = " + str(args.StartSweep) + "\n")
    MetaOutFD.write("#EndSweep = " + str(args.EndSweep) + "\n")
    MetaOutFD.write("#StartSamp = " + str(args.StartSamp) + "\n")
    MetaOutFD.write("#EndSamp = " + str(args.EndSamp) + "\n")
    MetaOutFD.write("#Scale = " + str(args.Scale) + "\n")
    MetaOutFD.write("#Log = TRUE\n")

## Prepare to loop

    StackCenter = int(numpy.fix((args.IncoDepth/2+0.5)))

    idxinco = 0
    for Stacked in read_and_stack_RADnh3(args.InputName,args.ChannelNum,args.SweepLength,args.StackDepth):
        truncStacked=Stacked[0:args.truncSweepLength]
        Dechirped = denoise_and_dechirp(truncStacked, Rchirp,args.blanking,args.truncSweepLength)
        # Get the magnitude
        if (idxinco == 0):
            IncoherentStack = numpy.abs(Dechirped)
            if args.PhsName:
                PhsStack = numpy.angle(Dechirped)
        else:
            IncoherentStack = numpy.vstack((IncoherentStack,numpy.abs(Dechirped)))
            if args.PhsName:
                PhsStack = numpy.vstack((PhsStack,numpy.angle(Dechirped)))

        idxinco = idxinco + 1
        if (idxinco == args.IncoDepth):
            idxinco = 0
            Incoherent = numpy.mean(IncoherentStack, axis=0)
            # Get just the center phase
            if args.PhsName:
                Phase = numpy.int32(PhsStack[StackCenter,...] * 16777216)
                Phase.byteswap(True)
                Phase.tofile(PhsOutFD)

            ScaledMag = numpy.int32(args.Scale * numpy.log10(Incoherent))
            ScaledMag.byteswap(True)
            ScaledMag.tofile(MagOutFD)
            if (enable_meta_index):
                MetaOutFD.write(str(record) + "\n")
            record=record+record_increment

if __name__ == "__main__":
    main(sys.argv)
