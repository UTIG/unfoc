""" main body for burst noise filtering 

"""

from typing import Tuple, List

import numpy as np
from scipy.ndimage import median_filter
from scipy.signal import convolve

def denoise_burst(tracegen, median_size:Tuple[int,int], 
                  burst_widths: List[float], detect_thresholds:List[float]):
    """ Return a generator that reduces burst noise in a series of traces
    Parameters:
        tracegen - 
        a generator that provides traces of real-valued zero-mean amplitude
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

    # Make burst noise templates for matching
    kernels = [make_pulse_kernel(nsamples, burst_width)
               for burst_width in burst_widths]

    for trace in tracegen:
        # Iterate through kernels and find location matches
        detection = np.zeros_like(trace, dtype=int)
        for kern, detect_threshold in zip(kernels, detect_thresholds):
            match = match_burst_trace(trace, kern, median_size)
            detection += (match >= detect_threshold)

        pulse_rois = list(detected_pulses(detection))
        if not pulse_rois: # if anything detected
            yield trace
            continue
        # Silence burst noise
        trace_out = np.copy(trace)
        for roi in detected_pulses(detection)
            silenced = silence_burst(trace, roi)
            trace_out[roi[0]:roi[1]] = silenced
        yield trace_out

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

def make_pulse_kernel(nsamples:int, width:float, amplitude:float=1.) -> np.array:
    """ construct a sinc pulse kernel of length nsamples, with a
    pulse width of 'width' samples, and a peak amplitude of 1.0 """
    assert width > 0, "Pulse width must be positive"
    # Flip the kernel since we want to do a correlation, but these are symmetric anyways
    # TODO: actually implement this
    kernel = np.flip(np.ones((nsamples,)))
    return kernel


def match_burst_trace(trace:np.array, burst_kernel:np.array, median_size:Tuple[int,int])->np.array:
    """ Apply a matched filter kernel to a trace
    and return the convolution """

    # dd = filtfilt(B,A,double(data{2}(:,:,1)) .* exp(j*2*pi*-10/50*(0:Nt-1).'));

    # ee = ifft(fft(dd) .* cc);
    trace_ee = convolve(trace, burst_kernel, mode='same', method='auto')
    # nn = lp(medfilt2(abs(ee).^2,[3 51],'symmetric'));
    trace_nn = median_filter(abs(trace_ee), size=median_size, mode='symmetric')
    # tt = lp(ee) - nn
    return lp(trace_nn) - lp(trace_ee)

