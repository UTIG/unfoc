unfoc Module Design
===================

2022-03-08 Gregory Ng

This document provides some discussion regarding how the unfoc module
is organized.

The previous iteration of the module simultaneously read all of
the input radar data, demultiplexed the radar channels, and wrote
filtered data to disk.

The old code was undoubtedly very difficult to understand and
kept track of a lot of state.

The current iteration processes one output data product at a time.
This has the potential to be somewhat slower, but because we are
doing less, the code is significantly simpler.

While it might seem slower and less efficient to read and process
one channel at a time, the theoretical time advantage/disadvantage is
unclear.  On one hand, you are rereading the same data more than once;
but on the other hand, disk caching minimize additional work done,
and parallel processing is easily achieved.

Because this processing appears somewhat I/O limited, there are no
plans to allow additional block-level parallelism (which seems to
make sense for 1D and especially 2D focusing).

primary entry point
-------------------

The primary entry point for the module is the `unfoc` function.

This function takes the usual processing parameters for inputs and
outputs, and then calls `unfoc_chan`, which runs the processing
code for each requested output data channel, assigning each channel
to among a pool of one or more workers.

processing pipeline
-------------------

`unfoc_chan` sets up a processing pipeline that creates a source
of raw radar traces.  It then chains this source into a coherent
stacker, a dechirp function, and a which is then passed through
functions for coherent stacking, dechirp, and incoherent stacking.

After the pipeline is set up, an output handler takes care of getting
each incoherent trace and formatting it for output, and writing to disk.

This pipeline is effectively driven by the output, and input traces
are pulled as needed to create the output incoherent traces, and
when there are no more available, it quits.

Input data readers
------------------

The input data reader classes provide an easy-to-use interface to
bxds files in python, and hide the task of parsing the bxds files.

Testing
-------

In the `regress` directory, run `run_coverage.sh` for regression
testing.  This script runs individual unit test python scripts.

References
----------

- Peters, M. E., Blankenship, D. D., & Morse, D. L. (2005).
  *Analysis techniques for coherent airborne radar sounding: Application to West Antarctic ice streams*.
  **Journal of Geophysical Research: Solid Earth**, 110(B6), B06303.
  https://doi.org/10.1029/2004JB003222
