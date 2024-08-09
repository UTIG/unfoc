""" main body for burst noise filtering 

"""

from typing import Tuple, List
import logging
from pathlib import Path

import numpy as np
from scipy.ndimage import median_filter
from scipy.signal import convolve, butter, filtfilt, medfilt2d, medfilt
#import matplotlib.pyplot as plt

from ..read import Trace

def denoise_burst(tracegen, median_size:Tuple[int,int], 
                  burst_widths: List[float], detect_thresholds:List[float]):
    """ Return a generator that reduces burst noise in a series of traces
    Parameters:
        tracegen - 
        a generator that provides Trace objects of real-valued zero-mean amplitude
        data (such as from a bxds file).  These can be 16-bit integer data

        median_size
        A tuple indicating the size of the 2D median filter for detecting background level.
        Dimensions are (slow time, fast time).  A reasonable value is (3, 51)

        Currently, only a median filter slow time dimension of 1 is supported
        (no forward or aft time dependence)
       
        burst_widths: List[float]
        A sequence of widths of bursts to use for the matched filter, in fast time samples.
        A reasonable value for this is [20], or perhaps [5, 15, 30]

        detect_threshold:
        A list of detection threshold values for burst detection in dB.
        The length of this sequence must be the same as for burst_widths.
        25 is a reasonable value

    Returns
        Yields a sequence of traces with the same length as the input sequence,
        where the outputs are filtered to remove the burst noise
    """

    assert median_size[0] == 1, "Median must be 1D for now"
    assert median_size[0] % 2 == 1 and median_size[1] % 2 == 1, \
           "Median filter dimensions must be odd"
    assert min(burst_widths) > 0, "Burst widths must be positive"
    assert len(burst_widths) == len(detect_thresholds), \
           "Number of burst widths and detection thresholds must match"

    if False:
        # Make burst noise templates for matching
        kernels = [make_pulse_kernel(burst_width, nsamples=None, amplitude=1000.)
                   for burst_width in burst_widths]
    else:
        kernels = [load_pulse_kernel()]
    numpulses = 0
    lpf = butter(2, 11/25) # downconversion lowpass filter
    expon = None # downconversion complex exponential
    logging.info("detect thresholds: %r", detect_thresholds)
    for ii, traceobj in enumerate(tracegen):
        trace = traceobj.data
        # Iterate through kernels and find location matches
        detection = np.zeros_like(trace, dtype=np.int8)

        # trace downconverted to baseband
        if expon is None:
            Nt = len(trace)
            # dd = filtfilt(B,A,double(data{2}(:,:,1)) .* exp(j*2*pi*-10/50*(0:Nt-1).'));
            expon = np.exp(1j*2*np.pi*-10/50*np.arange(Nt))
        trace_bb = filtfilt(*lpf, trace.astype(float) * expon)

        # TODO: refactor this to allow all of the median filters
        # to be done at once.  This should be more efficient since
        # the setup for the median_filter function only needs to be
        # called once for the entire 3-dimensional array
        for kern, detect_threshold in zip(kernels, detect_thresholds):
        
            match = match_burst_trace(trace_bb, kern, median_size)
            assert match.shape == detection.shape, "unexpected match shape"
            #detection += (match >= detect_threshold)
            np.maximum(detection, (match >= detect_threshold), out=detection)


        pulse_rois = list(detected_pulses(detection))
        if not pulse_rois: # if nothing detected
            yield traceobj
            continue
        numpulses += len(pulse_rois)
        # Silence burst noise
        trace_out = np.copy(trace)
        for roi in pulse_rois:
            trace_out[roi[0]:roi[1]] = silence_burst(trace, roi)
        yield Trace(traceobj.channel, trace_out, traceobj.ct)

        #logging.warning("trace %d Number of Pulse ROIs: %d", ii, numpulses)
        """ debugging plots
        if ii == 932:
            plt.clf()
            plt.plot(lp(trace), alpha=0.5, label='orig')
            plt.plot(lp(trace_out), alpha=0.5, label='silenced')
            plt.plot(lp(trace) - lp(trace_out), label='diff')
            plt.grid()
            plt.legend()
            plt.show()
        """

    #logging.debug("Number of Pulse ROIs: %d", numpulses)

def detected_pulses(detection:np.array, gap:int=14):
    """ Given a detection array, return a sequence of
    contiguous intervals

    detection is an array of integers where any nonzero
    value is considered a detection.

    gap is the minimum gap in detections to trigger
    a noncontiguous region.

    ROIs will be expanded by gap//2, up to the
    min/max dimensions of the detection array.
    TODO: use median filter and make this shorter/quicker
    """
    start_idx = None # first detection in this group
    prev_idx = None # previous detection index
    for ii, v in enumerate(detection):
        if v: # detected
            prev_idx = ii # regardless of being in detection, set prev
            if start_idx is None: # we are not in a detection region
                start_idx = ii
        else: # not detected
            if start_idx is not None and ii - prev_idx >= gap:
                yield expand_roi(start_idx, prev_idx+1, gap, detection)
                start_idx, prev_idx = None, None # Reset state

    if start_idx is not None: # Open interval hit the end, return
        yield expand_roi(start_idx, prev_idx+1, gap, detection)

def expand_roi(idx0:int, idx1:int, gap:int, detection:np.array)->Tuple[int,int]:
    """ idx0 is the start index,
    idx1 is the index of the first non-detection
    (i.e., a half-open interval)
    """
    roi_e_start = max(0, idx0-gap//2)
    roi_e_end = min(len(detection), idx1+gap//2)
    return roi_e_start, roi_e_end


def silence_burst(trace:np.array, roi:Tuple[int,int])->np.array:
    """ Silence the burst by taking the ends and filling them toward
    the middle
    TODO: make a separate function that does a linear interpolation
    rather than taking beginning and end
    """
    length = roi[1] - roi[0]
    arr = np.empty((length,), dtype=trace.dtype)
    mid = length // 2
    arr[0:mid] = trace[roi[0]]
    arr[mid:] = trace[roi[1]-1]
    return arr

def make_pulse_kernel(width:float,amplitude:float=1., nsamples:int=None) -> np.array:
    """ construct a sinc pulse kernel of length nsamples, with a
    pulse width of 'width' samples, and a peak amplitude of 1.0 
    Recommend making the kernel an odd number of samples for symmetry
    """
    assert width > 0, "Pulse width must be positive"
    if nsamples is None:
        nsamples = int(np.ceil(width*2.25))
        nsamples += ((nsamples+1) % 2) # enforce odd number of samples

    # Generate a sinc pulse with width 'width' samples in a kernel of nsamples samples
    w = nsamples / width
    x = np.linspace(-w, w, num=nsamples)
    kernel = np.sinc(x)*amplitude
    return kernel

def load_pulse_kernel():
    filename = Path(__file__).parent / 'burst_noise_PPT_MKB2o_Y01e_trace932_samp2985.npz'
    kern = np.load(filename)['burst_noise_sample']
    return kern#return np.flip(kern)





def match_burst_trace(trace_bb:np.array, burst_kernel:np.array, median_size:Tuple[int,int])->np.array:
    """ Apply a matched filter kernel to a trace
    and return the convolution 
    Assumes the trace has already been downconverted to baseband
    and that the kernel is also already adownconverted complex exponential
    """

    # ee = ifft(fft(dd) .* cc);
    trace_ee = np.atleast_2d(convolve(trace_bb, burst_kernel, mode='same', method='auto'))

    # medfilt2d is faster, around 290 seconds as median_filter for median_size=(1,51)
    # medfilt is the same as median_filter, 377 seconds
    # nn = lp(medfilt2(abs(ee).^2,[3 51],'symmetric'));
    #trace_nn = median_filter(abs(trace_ee), size=median_size, mode='reflect')
    trace_nn = medfilt2d(abs(trace_ee), kernel_size=list(median_size))
    #trace_nn = medfilt(abs(trace_ee[0,:]), kernel_size=median_size[1])
    # tt = lp(ee) - nn
    tt = lp(trace_ee/trace_nn)[0,:] #abs(lp(trace_nn) - lp(trace_ee))
    assert tt.shape == trace_bb.shape, "Input shape doesn't match output: input.shape=%r tt.shape=%r" % (trace.shape, tt.shape)
    return tt


def lp(trace:np.array):
    """ Log power of an amplitude trace """
    return 20*np.log10(abs(trace))
