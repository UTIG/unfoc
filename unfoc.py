#!/usr/bin/env python
#
# unfoc.py: Unfocused Radar Processor
# Output file for pik1 (4-byte signed integer, network order)
# Dechirp the JPL and MIRS radars.
#
# All output files are 4-byte network-order.
#
## (DONE): remove queueing because cache misses cause it to be slow.

import os
import gzip
import sys
import argparse
import numpy
import Queue
import threading
import logging
import struct
from collections import namedtuple
from scipy import signal
import cProfile

from parse_channels import parse_channels

################################################
# Enable metadata index writing.  
# This should normally be disabled only for legacy testing.
# It does not make things slower to output this data.
enable_meta_index=True
################################################

# from http://stackoverflow.com/questions/3160699/python-progress-bar
def update_progress(n,total):
    progress=float(n)/float(total)
    barLength = 20 # Modify this to change the length of the progress bar
    status = ""
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
        status = "error: progress var must be float\r\n"
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    if progress >= 1:
        progress = 1
        status = "Done...\r\n"
    block = int(round(barLength*progress))
    text = "\rPercent: [{0}] {1:.1f}% {2}/{3} {4}".format( "#"*block + "-"*(barLength-block), progress*100, n, total, status)
    sys.stdout.write(text)
    sys.stdout.flush()

class IncoStackState:
    """ This state is separated out from IncoStackFilter so we can use yield or as a functional call """
    def __init__(self, StackDepth, truncSweepLength, bDoPhs=False):
        self.StackDepth = StackDepth
        self.IncoherentStack = StackState1(StackDepth, truncSweepLength, False) # magnitude
        self.PhsStack        = StackState1(StackDepth, truncSweepLength, False) # phase
        self.bDoPhs = bDoPhs
        self.StackCenter = int(numpy.fix( (self.StackDepth/2+0.5) ))
        #print 'StackCenter {}'.format(self.StackCenter)
        
    def dostack(self, Dechirped):
        """ Trace is a complex, dechirped trace """
        mag = self.IncoherentStack.dostack( numpy.abs( Dechirped[1]), Dechirped[2])

        if self.bDoPhs:
            phs = self.PhsStack.dostack( numpy.angle( Dechirped[1] ), Dechirped[2])
        if mag != None:
            mag = numpy.mean( mag[0]  , axis=0 )
            if self.bDoPhs:
                phs = phs[0][self.StackCenter,...]
            else:
                phs = None
            return (Dechirped[0],mag, phs,Dechirped[2])

class QueueSource:
    """ Gets records from a python Queue to be used as an input to a filter """
    def __init__(self, maxsize=0):
        self.q = Queue.Queue(maxsize)
        pass

    def records(self):
        """ yield records.  If records is None, then quit. """
        while True:
            item = self.q.get()
            if item != None:
                yield item
                self.q.task_done()
            else:
                return

    def put(self, item, block=True, timeout=None):
        """ Add items to the queue for processing (to be returned by records() ) """ 
        return self.q.put(item, block, timeout)

def cinterp(Data, index):
    r = (numpy.abs(Data[index-1]) + numpy.abs(Data[index+1])) / 2
    t1 = numpy.angle(Data[index-1])
    t2 = numpy.angle(Data[index+1])
    if (numpy.abs(t1 - t2) > numpy.pi):
        t1 = t1 + 2 * numpy.pi
    theta = (t1 + t2) / 2
    Data[index] = r * (numpy.cos(theta) + 1j * numpy.sin(theta))
    return Data

def denoise_and_dechirp(Stacked, Rchirp, blanking, truncSweepLength, bDoCinterp = True):
    Stacked[0:blanking] = numpy.zeros(blanking)

#find peak energy below blanking samples
## [n,m]=sort(Stacked);
## shifter=abs((m(truncSweepLength)));
    shifter=int(numpy.median(numpy.argmax(Stacked)));
    Stacked=numpy.roll(Stacked,-shifter);
        
    DFT = numpy.fft.fft(signal.detrend(Stacked))

    if bDoCinterp:
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

def denoise_and_dechirp_gen(cohstacks, Rchirp, blanking, truncSweepLength, bDoCinterp = True):
    for trace in cohstacks:
        Dechirped = denoise_and_dechirp(trace[1][0:truncSweepLength], Rchirp, blanking, truncSweepLength, bDoCinterp)
        #print 'denoise_and_dechirped'
        #print trace
        yield (trace[0], Dechirped, trace[2])

class StackState1:
    """ State for stacking traces into blocks """
    def __init__(self, StackDepth, SweepLength, presum=True):
        self.StackDepth = StackDepth
        self.stacks = numpy.zeros((StackDepth, SweepLength), dtype=numpy.float64)
        self.idx = 0
        self.count = 0
        self.presum = presum

    def dostack(self,sweep,seq):
        self.stacks[self.idx,...] = sweep
        self.idx   += 1
        self.count += 1
        if self.idx >= self.StackDepth:
            self.idx = 0
            if self.presum:
         #       print 'StackState1.dostack'
         #       print (self.stacks.sum(axis=0),seq-(self.StackDepth/2))
                return (self.stacks.sum(axis=0),seq-(self.StackDepth/2))
            else:
                return (self.stacks,seq)
        else:
            return None

class StackState:
    """ Accepts multiple completed stacks and returns a stack if both input channel queues contain data """
    def __init__(self, ChannelSpec, StackDepth):
        self.ChannelSpec = ChannelSpec
        self.StackDepth  = StackDepth

        self.fullstacks0 = []
        self.fullstacks1 = []
        self.seq0 = []
        self.seq1 = []

    def getstack(self):
        # Logic lifted from read_and_stack.cc
        bDo0 = (self.ChannelSpec.chan0in > 0)
        bDo1 = (self.ChannelSpec.chan1in > 0)
        bReady0 = bDo0 and len(self.fullstacks0) > 0
        bReady1 = bDo1 and len(self.fullstacks1) > 0
        scale0 = self.ChannelSpec.scalef0 / float(self.StackDepth)
        scale1 = self.ChannelSpec.scalef1 / float(self.StackDepth)
        array_out1 = None
        seq_out=None

        if not ((bDo0 and not bReady0) or (bDo1 and not bReady1)) :
            # assert self.fullstacks0[-1].shape[0] == self.StackDepth?

            if bDo0 and bDo1:
                if not (self.seq0 == self.seq1):
                    logging.info('SEQ Mismatch: ch0={} ch1={}'.format(self.seq0,self.seq1))
                array_out1 = scale0 * self.fullstacks0[-1] + scale1 * self.fullstacks1[-1]
                seq_out=self.seq0[-1]
            elif bDo0 and not bDo1:
                array_out1 = scale0 * self.fullstacks0[-1]
                seq_out=self.seq0[-1]
            elif not bDo0 and bDo1:
                array_out1 = scale1 * self.fullstacks1[-1]
                seq_out=self.seq1[-1]

            if bReady0: 
                self.fullstacks0.pop();
                self.seq0.pop();
            if bReady1: 
                self.fullstacks1.pop();
                self.seq1.pop();
            return (self.ChannelSpec.chanout,array_out1,seq_out)
 
    def addstack( self, chan, stack ):
        if len(self.fullstacks0) >= 1000 or len(self.fullstacks1) >= 1000:
            raise Exception("overflow: p1cs={4:s} fullstacks0={0:d} fullstacks1={1:d} count0={2:d} count1={3:d}".format(len(self.fullstacks0), len(self.fullstacks1), self.count0, self.count1, self.ChannelSpec) )
        if self.ChannelSpec.chan0in == chan:
            self.fullstacks0.insert(0,stack[0])
            self.seq0.insert(0,stack[1])
        elif self.ChannelSpec.chan1in == chan:
            self.fullstacks1.insert(0,stack[0])
            self.seq1.insert(0,stack[1])
        else:
            # no change in state, so no need for getstacks
            return None

        return self.getstack()

# Yields a tuple of an output channel and a NDArray of values
# coherently stacked
# tracegen should be a sequence (generator) of tuples where traces[0] is the channel id, 
# traces[1] is an ndarray that is the trace
# and traces[2] is the sequnece number of the trace
def Stacks_gen(traces, ChannelSpecs, StackDepth, SweepLength=3437):
    stackstates = []
    ss0 = {}
    for (i,p1cs) in enumerate(ChannelSpecs):
        logging.debug("p1cs[%d]=%s"% (i,str(p1cs)))
        stackstates.append(StackState(p1cs, StackDepth))
        for chan in (p1cs.chan0in, p1cs.chan1in):
            if chan > 0 and chan not in ss0:
                ss0[chan] =  StackState1(StackDepth, SweepLength)

    ## TODO: make this read function selectable to accommodate for RADnh4
    #for trace in read_RADnh3_gen(self.InputName, self.ChannelSpecs, self.SweepLength):
    for trace in traces:
        # Stack into appropriate bins
        # Iterate over all the channel specs and see which ones need to be stacked

        if trace[0] in ss0:
            cohstack = ss0[ trace[0] ].dostack(numpy.float64(trace[1]),trace[2])

            if cohstack != None:
                for (i,p1cs) in enumerate(ChannelSpecs):
                    # Test to prevent extraneous calls to addstack()
                    if p1cs.chan0in == trace[0] or p1cs.chan1in == trace[0]:
                        result = stackstates[i].addstack(trace[0],cohstack)
                        # If this stack state yielded a full stack, yield it
                        if result != None:
                            yield result



# Yields a tuple of an output channel and a NDArray of values
# Incoherently stacked
# tracegen should be a sequence (generator) of tuples where traces[0] is the channel id, 
# and traces[1] is an ndarray that is the trace
def IncoStacks_gen(traces, ChannelSpecs, StackDepth, truncSweepLength=3200, bDoPhs=False):
    ss0 = {}
    for (i,p1cs) in enumerate(ChannelSpecs):
        chan = p1cs.chanout
        if chan > 0 and chan not in ss0:
            ss0[chan] = IncoStackState(StackDepth, truncSweepLength, bDoPhs)

    ## TODO: make this read function selectable to accommodate for RADnh4
    #for trace in read_RADnh3_gen(self.InputName, self.ChannelSpecs, self.SweepLength):
    for trace in traces:
        # Stack into appropriate bins
        # Iterate over all the channel specs and see which ones need to be stacked

        if trace[0] in ss0:
            stack = ss0[ trace[0] ].dostack( trace )

            if stack != None:
                yield stack
                            
# Read individual traces out of RADnh3 file        
def read_RADnh3_gen(InputName, ChannelSpecs, ct, SweepLength=3437):
    radnh3_header_t = namedtuple('radnh3_header', 'nsamp nchan vr0 vr1 choff resvd1 resvd2')
    n=0

    # Construct a list of all channel offsets we want
    choffs = set()
    for p1cs in ChannelSpecs:
        choffs.add(((p1cs.chan0in-1)//2)*2)
        choffs.add(((p1cs.chan1in-1)//2)*2)


    with open(InputName, 'r') as fd:
        while True:
            # read header
            buff = fd.read(8)
            if len(buff) < 8:
                return
            radnh3_header = radnh3_header_t._make(struct.unpack_from('>HBBBBBB', buff ))

            # Legacy value of 0xff is equivalent to offset 0
            if radnh3_header.choff == 0xff:
                choff = 0
            else:
                choff = 0x0f & radnh3_header.choff

            if choff in choffs:
                traces = numpy.fromfile(fd, dtype='>i2', count=SweepLength*2)
                traces.shape = (2, SweepLength)
                yield (choff+1, traces[0][...],ct['seq'][n])
                yield (choff+2, traces[1][...],ct['seq'][n])
            else:
                # Skip the rest of this record without yielding any values
                fd.seek(4*SweepLength, os.SEEK_CUR)
            if (n % 2400) == 0: 
                update_progress(n,ct.shape[0])
            n += 1

# Read individual traces out of RADjh1 files        
def read_RADjh1_gen(InputName, ChannelSpecs, ct, SweepLength=3200):
    n=0

    # Construct a list of all channel offsets we want
    with open('{}1'.format(InputName), 'r') as fd1, open('{}2'.format(InputName), 'r') as fd2:
        ### Redundant.  GNG
        while True and n < ct.shape[0]:
            traces1 = numpy.fromfile(fd1, dtype='<i2', count=SweepLength)
            traces2 = numpy.fromfile(fd2, dtype='<i2', count=SweepLength)
            yield (1, traces1 ,ct['seq'][n])
            yield (2, traces2 ,ct['seq'][n])
            if (n % 2400) == 0: 
                update_progress(n,ct.shape[0])
            n=n+1


class DataFilterBase:
    def __init__(self, src):
        self.src = src

    def records(self):
        pass


class ComplexOutputFile(DataFilterBase):
    pass

class  PIK1OutputFile(ComplexOutputFile):
    """ Outputs magnitude and phase information for NDArray that comes through.
        Optionally computes the mean before writing
        It can be used as a tee filter and tap output at any phase in the 
        processing stream.
    """
    # File descriptors
    MagOutFD  = None
    PhsOutFD  = None
    MetaOutFD = None
    
    # File names
    MagFileName  = None
    PhsFileName  = None
    MetaFileName = None

    IncoDepth = None
    MagScale  = None
    ChannelNum = None

    def __init__(self, src, MagScale):
        self.MagScale = MagScale
        ComplexOutputFile.__init__(self,src)

    def open(self,basepath, channel, InputFileName, args, bDoPhs=False, bDoIdx=True):
        self.close()

        # TODO: remove references to args as possible
        self.enable_meta_idx = bDoIdx
        self.ChannelNum = channel
        self.StackDepth = args.StackDepth
        self.IncoDepth = args.IncoDepth

        self.MagFileName = "{0:s}/{1:s}{2:d}".format(basepath,args.MagName, channel)
        self.PhsFileName = "{0:s}/{1:s}{2:d}".format(basepath, args.PhsName, channel)
        self.MetaFileName = "{0:s}/{1:s}{2:d}.meta".format(basepath, args.MetaName, channel)

        ## Open Input and Output files

        self.MetaOutFD = open(self.MetaFileName, 'w')
        self.MetaOutFD.write('#InputName = "' + InputFileName + '"\n')
        logging.info("writing %s" % self.MagFileName)
        self.MagOutFD = open(self.MagFileName, 'wb')
        self.MetaOutFD.write('#MagName = "' + self.MagFileName + '"\n')

        if bDoPhs:
            self.PhsOutFD = open(self.PhsFileName, 'w')
            self.MetaOutFD.write('#PhsName = "' + self.PhsFileName + '"\n')

        ###### FIXME? - request hdr file from xlob and read/forward
        # TODO: collect this information into a data structure rather than reading directly from args.
        ###### Make this Meta similar
        self.MetaOutFD.write("#ChannelNum = " + str(self.ChannelNum) + "\n")
        self.MetaOutFD.write("#StackDepth = " + str(self.StackDepth) + "\n")
        self.MetaOutFD.write("#IncoDepth = "  + str(self.IncoDepth) + "\n")
        self.MetaOutFD.write("#CenterMult = " + str(args.CenterMult) + "\n")
        self.MetaOutFD.write("#MaxDepth = "   + str(args.MaxDepth) + "\n")
        self.MetaOutFD.write("#StartSweep = " + str(args.StartSweep) + "\n")
        self.MetaOutFD.write("#EndSweep = "   + str(args.EndSweep) + "\n")
        self.MetaOutFD.write("#StartSamp = "  + str(args.StartSamp) + "\n")
        self.MetaOutFD.write("#EndSamp = "    + str(args.EndSamp) + "\n")
        self.MetaOutFD.write("#Scale = "      + str(self.MagScale) + "\n")
        self.MetaOutFD.write("#Log = TRUE\n")


        #self.record_increment=self.StackDepth*self.IncoDepth
        #self.record=self.record_increment/2

        # Write sys.ccn file, then close it.  
        StreamFileName='{0:s}/sys.ccn'.format(basepath)
        logging.debug('writing sys file: {0:s}'.format(StreamFileName))
        with open(StreamFileName,'w') as StreamFD:
            StreamFD.write(args.StreamName + "\n")
        

    def write_record(self, IncoStacked):
        # Write component files if enabled
        if self.PhsOutFD != None:
            Phase = numpy.int32(IncoStacked[2] * 16777216 )
            Phase.byteswap(True)
            Phase.tofile(self.PhsOutFD)

        if self.MagOutFD != None:
            ScaledMag = numpy.int32(self.MagScale * numpy.log10(IncoStacked[1]))
            ScaledMag.byteswap(True)
            ScaledMag.tofile(self.MagOutFD)

        if (self.enable_meta_idx):
            self.MetaOutFD.write(str(IncoStacked[3]) + "\n")
        #self.record += self.record_increment

    def records(self):
        bDoPhs = self.PhsOutFD != None
        for IncoStacked in self.src.records():
            self.write_record(IncoStacked)
            yield IncoStacked

    def close(self):
        # flush stacks,
        # warn about incomplete stacks
        # close file handles

        if self.MagOutFD != None:
            self.MagOutFD.close()
            self.MagOutFD = None

        if self.PhsOutFD != None:
            self.PhsOutFD.close()
            self.PhsOutFD = None

        if self.MetaOutFD != None:
            self.MetaOutFD.close()
            self.MetaOutFD = None


def worker_writer(p1file):
    logging.info("Starting channel %d" % p1file.ChannelNum)
    for rec in p1file.records():
        pass
    logging.info("Finished")


def threads_all_alive(list_threads):
    for t in list_threads:
        if not t.is_alive():
            return False
    return True



def main(argv):
    parser = argparse.ArgumentParser(description='Pulse compress radar data')
# filename of 2-byte radar file
    parser.add_argument('--InputName', required=True)
# filename of counter-timer file
    parser.add_argument('--InputCT', required=True)
# Channel Number to compress (counts from 1)
    parser.add_argument('--ChannelNum', default=0)
# Length of each sweep (in samples)
# TODO: This argument will soon be deprecated
    parser.add_argument('--SweepLength', required=True, type=int)
# Length of each output sweep (in samples)
    parser.add_argument('--truncSweepLength', required=True, type=int)

    parser.add_argument('--basepath', required=True, help='Base path for output files')

# GNG: these options are not used in unfoc
## filename of real part of raw compressed waveform
##    parser.add_argument('--IName')
## filename of imaginary part of raw compressed waveform
##    parser.add_argument('--QName')


# filename of rectified output
    parser.add_argument('--MagName')
# filename of phase output
    parser.add_argument('--PhsName')
    parser.add_argument('--MetaName',help='filename for meta data output')
    parser.add_argument('--StackDepth', required=True, type=int,help='coherent stacking depth for this output')
    parser.add_argument('--IncoDepth', required=True, type=int,help='incoherent stacking depth for this output')
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
    parser.add_argument('--Scale', type=int, default=20000, help='output scale factor: default is 1000*dB')
    parser.add_argument('--blanking', type=int, default=50, help='blank out (zero) this many samples at the top (beginning) of the record')
# Stream name.  Required to know data format
    parser.add_argument('--StreamName', required=True,help='Stream name for input data')
# bandpass sampling, false is for legacy hicars. true is for marfa
# Disables cinterp and flips the chirp
    parser.add_argument('--BandpassSamplingMode', action='store_true')
    parser.add_argument('--debug', action='store_true',help='Print debugging messages')
    parser.add_argument('--progress', action='store_true',help='Show progress bar')


    args = parser.parse_args()

    if args.debug:
        LOGLEVEL=logging.DEBUG
    else:
        LOGLEVEL=logging.INFO
    logging.basicConfig(level=LOGLEVEL,
                    format='{0:s}: %(relativeCreated)8d [%(levelname)-5s] (%(process)d %(threadName)-10s) %(message)s'.format(__file__),
                    )

    if args.EndSamp == 0 or args.EndSamp > args.SweepLength:
        args.EndSamp = args.SweepLength
    if args.StartSamp > args.EndSamp:
        sys.exit('{}: bad start and end samples.'.format(__file__))

    if args.StackDepth <= 0 or args.CenterMult <= 0 or args.MaxDepth <= 0:
        sys.exit('{}: all depths must be positive.'.format(__file___))

    if args.SweepLength <= 0:
        sys.exit('{}: sweep length must be positive'.format(__file__))

# Obtain reference chirp
    I = numpy.array([-63, -92, -109, -75, -87, -50, -116, -154, -22, 68, 141, -610, 1461, 3807, -6147, -5375, 10651, -4412, -9810, 15386, -3070, -14499, 15130, 3677, -15935, 3743, 13362, -5884, -13301, 8455, 12542, -8744, -11977, 5105, 13754, -961, -14342, -5184, 10294, 12194, -2709, -14352, -8807, 5965, 15350, 8368, -6605, -14990, -11515, 196, 11276, 15490, 10300, -645, -10730, -15307, -13379, -7342, 377, 7264, 11662, 13435, 11530, 3243, -4865, -4427, -3233, -4000, -2472, -2498, -2361, -1230, -1311, -618, -578, -569, -121, -319, 206, 328, 436, 613, 318, 514, 353, 277, 221, 34, 250, 132, 199, 189, 75, 190, 65, 106, 19, -64, -14, -117])
    if not args.BandpassSamplingMode:
        rchirp = numpy.flipud(I)
    else:
        rchirp = I
    Rchirp = numpy.fft.fft(rchirp,n=args.truncSweepLength)

# Compute Hamming Filter
    From = round(2.5 * args.truncSweepLength / 50)
    To = round(17.5 * args.truncSweepLength / 50)
    hamming = numpy.sin(numpy.linspace(0,1,num=To-From+1) * numpy.pi)
    hfilter = numpy.hstack((numpy.zeros(From), hamming*2, numpy.zeros(args.truncSweepLength-2*To-1), numpy.zeros(hamming.size), numpy.zeros(From-1)))

## Filter Chirp
    # disable for production work
    # plt.plot(hfilter, range(0,filter.size))
    Rchirp = numpy.multiply(Rchirp, hfilter)

    ChannelSpecs = parse_channels( args.ChannelNum )
    logging.debug(ChannelSpecs)

    # Read CT file
    logging.debug('reading ct file {0:s}'.format(args.InputCT))
    try:
        ct=numpy.loadtxt(args.InputCT,dtype={
            'names':['P','S','T','seq','YY','MM','DD','hh','mm','ss','fs','ct'],
            'formats':['|S8','|S8','|S8','int32','int16','int8','int8','int8','int8','int8','int8','int32']})
    except:
        logging.error('failed to read ct file')
        exit()
    logging.debug('processing radar data with {0:d} raw traces'.format(ct.shape[0]))

    # Read traces from file
    # TODO: this if structure could be made more inherity GNG
    if args.StreamName == 'RADnh3':
        tracegen = read_RADnh3_gen(args.InputName, ChannelSpecs, ct, args.SweepLength) 
    elif args.StreamName == 'RADjh1':
        tracegen = read_RADjh1_gen(args.InputName, ChannelSpecs, ct, args.SweepLength) 
    else:
        logging.error('Unknown stream name {0:s}'.format(args.StreamName))
        exit()
    

    #i=1
    #for trace in tracegen:
    #    if i < 10:
    #        print trace
    #        i=i+1
    #    else:
    #        exit()
    # Demultiplex stacks and generate coherently-stacked traces
    stackgen = Stacks_gen(tracegen, ChannelSpecs,args.StackDepth, args.SweepLength)
    # Dechirp coherent stacks
    dechirpgen = denoise_and_dechirp_gen(stackgen, Rchirp, args.blanking, args.truncSweepLength, not args.BandpassSamplingMode)

    # Incoherently stack 
    istackgen = IncoStacks_gen(dechirpgen, ChannelSpecs, args.IncoDepth, args.truncSweepLength, bDoPhs=args.PhsName != None)

    outfiles = {}
    #for p1cs in ChannelSpecs:
    #    outfiles[ p1cs.chanout ] = PIK1OutputFile( p1src, args.Scale )
#    outfiles[ p1cs.chanout ].open(args.basepath, p1cs.chanout, args.InputName, args, args.PhsName != "")

    bMultithread = False
    if bMultithread:
        # Start the writer threads.
        threads = []
        queues  = []

        for (i, p1cs) in enumerate(ChannelSpecs):
            qs1 = QueueSource(10)

            outfiles[ p1cs.chanout ] = PIK1OutputFile( qs1, args.Scale )


            outfiles[ p1cs.chanout ].open(args.basepath, p1cs.chanout, args.InputName, args, args.PhsName != None)

            # This thread pulls records from its source
            t = threading.Thread(target=worker_writer, args=(outfiles[ p1cs.chanout ],), name="ch%d" % ChannelSpecs[i].chanout)
            threads.append(t)
            queues.append(qs1)
            t.start()    
            try:
                # pump records from the source file into the queue
                for rec in istackgen:
                    try:
                        for (i, qs1) in enumerate(queues):
                            if rec[0] == ChannelSpecs[i].chanout:
                                qs1.put(rec, True, 1.0)

                    except Queue.Full:
                        if not threads_all_alive(threads):
                            break
                logging.info("Done writing.")
            finally:
                for q in queues:
                    try:
                        q.put(None)
                    except Queue.Full:
                        pass

        logging.debug("Waiting for threads to finish.")
        # wait for all queues to finish.
        for t in threads:
            t.join()

        logging.debug("Threads finished.")
        
    else:
        # Initialize output files
        for p1cs in ChannelSpecs:
            outfiles[ p1cs.chanout ] = PIK1OutputFile( None, args.Scale )
            outfiles[ p1cs.chanout ].open(args.basepath, p1cs.chanout, args.InputName, args, args.PhsName != None)
        for rec in istackgen:
            if rec[0] in outfiles:
                outfiles[ rec[0] ].write_record( rec)
            else:
                raise Exception("Invalid output channel %d" % rec[0] )

    update_progress(ct.shape[0],ct.shape[0])

    for p1cs in ChannelSpecs:
        outfiles[ p1cs.chanout ].close()


if __name__ == "__main__":
    bDoProfile = False
    if bDoProfile:
        import os
        prof_file = "/tmp/pyk1b.{0:d}.prof".format(os.getpid())
        cProfile.run('main(sys.argv)', prof_file)
        import pstats
        p = pstats.Stats(prof_file)
        p.sort_stats('cumulative').print_stats(50)
    else:
        main(sys.argv)
