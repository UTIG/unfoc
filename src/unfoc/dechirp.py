#!/usr/bin/env python3

"""
Math related to dechirping, local oscillator noise suppression and frequency domain filtering.

This module provides tools for:
- Removing coherent LO components from FFT spectra
- Dechirping traces using a reference chirp
- Generating reference chirps and frequency-domain filters
"""

import logging

import numpy as np
from scipy import signal

try:
    import typing
    from typing import Any, BinaryIO, Dict, Generator, List, Optional, Set, Tuple
except ImportError: #pragma: no cover
    pass  # Not installed on melt ...


def cinterp(sweep_fft, index):
    """ 
    Suppress local oscillator (LO) spike in an FFT spectrum by interpolation.

    Modifies sweep_fft in place to filter out coherent local oscillator
    component by interpolating in the complex fourier domain
    sweep_fft is fft of sweep, and index is a bin affected by the LO noise.

    Parameters
    ----------
    sweep_fft : np.ndarray
        Complex FFT array of the input radar sweep.
    index : int
        Index of the LO-contaminated bin to be replaced.
    """
    # type: (np.ndarray, int) -> None
    r = (np.abs(sweep_fft[index-1]) + np.abs(sweep_fft[index+1])) / 2
    t1 = np.angle(sweep_fft[index-1])
    t2 = np.angle(sweep_fft[index+1])
    if (np.abs(t1 - t2) > np.pi):
        t1 = t1 + 2 * np.pi
    theta = (t1 + t2) / 2
    # This is equivalent to r*np.exp(1j*theta), but it ain't broke.
    sweep_fft[index] = r * (np.cos(theta) + 1j * np.sin(theta))


# QUESTION: with numpy, is this modifying the input trace, or just the output?
def denoise_and_dechirp(trace, # type: np.ndarray
                        ref_chirp, # type: np.ndarray
                        blanking, # type: int
                        output_samples, # type: int
                        do_cinterp, # type: bool
                        detrend_type='linear'
                       ):
    """
    Denoise and dechirp a radar trace in the frequency domain.

    Applies blanking, time shift, detrending, FFT, optional LO suppression,
    and multiplies by a reference chirp in frequency domain. The result
    is inverse-transformed and realigned to yield the dechirped trace.

    Parameters
    ----------
    trace : np.ndarray
        Raw input radar trace to be processed.
    ref_chirp : np.ndarray
        Frequency-domain reference chirp used for dechirping.
    blanking : int
        Number of samples to zero out at the start (positive) or end (negative).
    output_samples : int
        Length of FFT and output trace.
    do_cinterp : bool
        Whether to apply LO spike suppression via `cinterp`.
    detrend_type : str, optional
        Detrending method passed to `scipy.signal.detrend`. Default is 'linear'.

    Returns
    -------
    np.ndarray
        The dechirped radar trace (complex-valued).
    """
    # type: (...) -> np.ndarray

    # Input trace is a 1 x output_samples array.
    if (blanking >= 0):
        trace[:blanking] = np.zeros(blanking)
    else:
        trace[-blanking:] = np.zeros(len(trace)+blanking)

    #find peak energy below blanking samples
    ## [n,m]=sort(trace);
    ## shifter=abs((m(output_samples)));
    # GNG -- this line wasn't doing what you think it was doing.
    #shifter = int(np.median(np.argmax(trace)));
    # You have to do this:
    #shifter = np.median(np.argwhere(listy == np.amax(listy)))
    shifter = int(np.argmax(trace))
    trace = np.roll(trace, -shifter);

    DFT = np.fft.fft(signal.detrend(trace, type=detrend_type))

    if do_cinterp:
        # Remove five samples per cycle problem
        n1 = int(np.round(output_samples * (1.0/5)))
        n2 = output_samples - n1
        cinterp(DFT, n1)
        cinterp(DFT, n2)
        # Remove the first harmonic for five samples
        n1 = int(np.round(output_samples * (2.0/5)))
        n2 = output_samples - n1
        cinterp(DFT, n1)
        cinterp(DFT, n2)

    # Do the dechirp
    Product = np.multiply(ref_chirp, DFT)
    Dechirped = np.fft.ifft(Product)
    Dechirped = np.roll(Dechirped, shifter)
    return Dechirped


def get_ref_chirp(bandpass, trunc_sweep_length):
    """
    Generate a reference chirp for frequency-domain dechirping.

    Parameters
    ----------
    bandpass : bool
        If True, return chirp as-is. If False, reverse chirp for matched filtering.
    trunc_sweep_length : int
        Desired length of the output FFT chirp.

    Returns
    -------
    np.ndarray
        Frequency-domain reference chirp.
    """

    # type: (bool, bool) -> np.ndarray
    I = np.array([-63, -92, -109, -75, -87, -50, -116, -154, -22, 68, 141,
                  -610, 1461, 3807, -6147, -5375, 10651, -4412, -9810, 15386,
                  -3070, -14499, 15130, 3677, -15935, 3743, 13362, -5884,
                  -13301, 8455, 12542, -8744, -11977, 5105, 13754, -961,
                  -14342, -5184, 10294, 12194, -2709, -14352, -8807, 5965,
                  15350, 8368, -6605, -14990, -11515, 196, 11276, 15490,
                  10300, -645, -10730, -15307, -13379, -7342, 377, 7264,
                  11662, 13435, 11530, 3243, -4865, -4427, -3233, -4000,
                  -2472, -2498, -2361, -1230, -1311, -618, -578, -569, -121,
                  -319, 206, 328, 436, 613, 318, 514, 353, 277, 221, 34, 250,
                  132, 199, 189, 75, 190, 65, 106, 19, -64, -14, -117])
    if not bandpass:
        rchirp = np.flipud(I)
    else:
        rchirp = I

    return np.fft.fft(rchirp, n=trunc_sweep_length)


def get_hfilter(trunc_sweep_length):
    """
    Generate a half-band sine window filter centered in a fixed band.

    The filter attenuates frequencies outside the 2.5–17.5 MHz range,
    scaled to the FFT length.

    Parameters
    ----------
    trunc_sweep_length : int
        FFT length, used to compute the sample indices corresponding to frequency cutoffs.

    Returns
    -------
    np.ndarray
        Frequency-domain filter array.
    """
    # type: (int) -> np.ndarray
    # Convert MHz to samples
    min_freq = int(round(2.5 * trunc_sweep_length / 50))
    max_freq = int(round(17.5 * trunc_sweep_length / 50))
    dfreq = max_freq - min_freq + 1
    hamming = np.sin(np.linspace(0, 1, num=dfreq) * np.pi)
    hfilter = np.hstack((np.zeros(min_freq),
                         hamming*2,
                         np.zeros(trunc_sweep_length - 2*max_freq - 1),
                         np.zeros(hamming.size),
                         np.zeros(min_freq - 1)))
    return hfilter
