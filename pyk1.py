#env python
#
# Dechirp the JPL and MIRS radars.
#
# All output files are 4-byte network-order.
#

import sys
import argparse
import struct

def main(argv):
    parser = argparse.ArgumentParser(description='Pulse compress radar data')
# filename of 2-byte radar file
    parser.add_argument('--InputName', required=true);
# Channel Number to compress (counts from 1)
    parser.add_arguemnt('--ChannelNum', default=0, type=int);
# Length of each sweep (in samples)
    parser.add_arguemnt('--SweepLength', required=true, type=int);
# Length of each output sweep (in samples)
    parser.add_arguemnt('--trunSweepLength', required=true, type=int);
# filename of real part of raw compressed waveform
    parser.add_arguemnt('--IName');
# filename of imaginary part of raw compressed waveform
    parser.add_arguemnt('--QName');
# filename of recitified output
    parser.add_arguemnt('--MagName');
# filename of phase output
    parser.add_arguemnt('--PhsName');
# filename for meta data output
    parser.add_arguemnt('--MetaName');
# coherent stacking depth for this output
    parser.add_arguemnt('--StackDepth', required=true, type=int);
# incoherent stacking depth for this output
    parser.add_arguemnt('--IncoDepth', required=true, type=int);
# stacks are centered on (or before) multiples of CenterMult
    parser.add_arguemnt('--CenterMult', required=true, type=int);
# Only output a sweep if MaxDepth sweeps could be stacked.
    parser.add_arguemnt('--MaxDepth', required=true, type=int);
# Output only sweeps with centers on or after this sweep
    parser.add_arguemnt('--StartSweep', type=int);
# Output only sweeps with centers before or on this sweep
    parser.add_arguemnt('--EndSweep', type=int);
# Output sweep samples starting with this one
    parser.add_arguemnt('--StartSamp', type=int);
# Output sweep samples ending with this one
    parser.add_arguemnt('--EndSamp', type=int);
# Output scale default is 1000*dB
    parser.add_arguemnt('--Scale', type=int, default=20000);
# Samples at the top of the record to blank out
    parser.add_arguemnt('--blanking', type=int, default=50);
# Stream name
    parser.add_arguemnt('--StreamName', required=true);

    if StartSamp < 1:
        StartSamp = 1
    if EndSamp > SweepLength:
        EndSamp = SweepLength
    if StartSamp > EndSamp:
        sys.exit('pyk1: bad start and end samples.')

    if StackDepth <= 0 or CenterMult <= 0 or MaxDepth <= 0
        sys.exit('pyk1: all depths must be positive')

    if SweepLength <= 0
        sys.exit('pyk1: sweep length must be positive');

    record_increment=StackDepth*IncoDepth
    record=round(record_increment/2)

# Obtain reference chirp
    I = numpy.fromFile('I.bin', dtype='>i4');
    Q = numpy.fromFile('Q.bin', dtype='>i4');

    chirp = complex(I, Q)
    rchirp = flipud(chirp)
    Rchirp = fft(rchirp,n=truncSweepLength)

# Compute Hamming Filter
    filter = zeros((truncSweepLength,1))
    from = round(2.5 * truncSweepLength / 50)
    to = round(17.5 * truncSweepLength / 50)
    diff = to - from
    hamming = sin(arange(0:1+1/diff:1/diff) * np.pi)  # arange is open on right
    filter(from:to) = hamming # illegal in python
    filter(truncSweepLength-to+2:truncSweepLength-from+2) = hamming # illegal in python

## Filter Chirp

    Rchirp = Rchirp dot filter

## Open Input and Output files

    MetaOutFD = open(MetaName)
    MetaOutFD.write('#MagName = "' + InputName + '"\n')

    MagOutFD = open(MagName)
    MetaOutFD.write('#MagName = "' + MagName + '"\n')

    PhsOutFD = open(PhsName)
    MetaOutFD.write('#PhsName = "' + PhsName + '"\n')

    ###### FIXME - request hdr file from xlob and read/forward
    ###### Make this Meta similar
    MetaOutFD.write("#ChannelNum = " + str(ChannelNum) + "\n");
    MetaOutFD.write("#StackDepth = " + str(StackDepth) + "\n");
    MetaOutFD.write("#IncoDepth = " + str(IncoDepth) + "\n");
    MetaOutFD.write("#CenterMult = " + str(CenterMult) + "\n");
    MetaOutFD.write("#MaxDepth = " + str(MaxDepth) + "\n");
    MetaOutFD.write("#StartSweep = " + StartSweep) + "\n");
    MetaOutFD.write("#EndSweep = " + str(EndSweep) + "\n");
    MetaOutFD.write("#StartSamp = " + StartSamp) + "\n");
    MetaOutFD.write("#EndSamp = " + str(EndSamp) + "\n");
    MetaOutFD.write("#Scale = " + str(Scale) + "\n");
    MetaOutFD.write("#Log = TRUE\n");

## Prepare to loop

    Done = False
    StackCenter = fix((StackDepth/2+0.5))
    Stack = zeros(truncSweepLength,StackDepth);
    end_of_data = -9999999.0

    while not Done:
        for idxinco in range(0,IncoDepth)
        for idxinco = 1:IncoDepth

    Stacked=read_and_stack_RADnh3(InputName,ChannelNum,SweepLength,StackDepth,end_of_data);

	if (Stacked(1) == end_of_data)
            fprintf(stderr, "dechirp.m: no full stacks left\n");
            Done = 1;
            break;
        else
            truncStacked=Stacked(1:3200);
            Dechirped = denoise_and_dechirp(truncStacked, Rchirp,blanking,truncSweepLength);
            % Get the magnitude
            IncoherentStack(1:truncSweepLength,idxinco) = abs(Dechirped);
            PhsStack(1:truncSweepLength,idxinco) = angle(Dechirped);
        endif	

    endfor

    if (!Done)
        % GNG: Why add 50/2?
        %Incoherent=mean(IncoherentStack+(50/2),2);
        Incoherent = mean(IncoherentStack, 2);
        % Get just the center phase
        Phase = PhsStack(1:truncSweepLength, StackCenter);
        LogMag = log10(Incoherent);
        ScaledMag = Scale * LogMag;
        MyWrite(MagOutFD, ScaledMag(StartSamp:truncSweepLength));
        MyWrite(PhsOutFD, Phase(StartSamp:truncSweepLength)*16777216);
   % 	fprintf(MetaOutFD,"%d\n", record);
	record=record+record_increment;
    endif

endwhile


if (MagOutFD != -1) fclose(MagOutFD); fclose(PhsOutFD); endif;
retval = 1;
endfunction

################################################################################
################################################################################
################################################################################

if __name__ == "__main__":
    main(sys.argv)
