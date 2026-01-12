# Unfocused Radar Processor

This package implements unfocused radar processing for radar data
collected by the UTIG/SOAR family of ice penetrating radar systems.

The primary purpose is to read radar data from disk and output an
unfocused data product (aka pik1).

It includes a module for reading radar data from disk and loading
into numpy as an array with timing metadata.


## More Detailed Documentation (Sphinx sources)

- [Usage Guide](docs/source/usage_guide.rst)
- [Installation](docs/source/installation.rst)
- [Package design](docs/source/package_design.rst)
- [Modules](docs/source/modules.rst)



# Command Line Usage

You can invoke the module as a script to process input data on disk
and produce the output unfocused data to another folder.

```
run_unfoc.py -i $WAIS/orig/xlob/ICP3/JKB2d/F56T01c/RADnh3/bxds -o /tmp/outdir \
             --stackdepth 10 --incodepth 5 --channels LoResInco1,LoResInco2
```

This will read the bxds file at `$WAIS/orig/xlob/ICP3/JKB2d/F56T01c/RADnh3/bxds`
(and an associated `ct` or `ct.gz` metadata file alongside it, and output
the requested output channels *LoResInco1 and *LoResInco2

To run each output channel on a separate core, add the `-j` or `--jobs` option:

```
run_unfoc.py -i $WAIS/orig/xlob/ICP3/JKB2d/F56T01c/RADnh3/bxds -o /tmp/outdir \
             --stackdepth 10 --incodepth 5 --channels LoResInco1,LoResInco2 \
             -j 2
```



# Python Usage

This same functionality can be called natively from python in largely the same way.

```python
from unfoc import unfoc

channels = 'LoResInco1,LoResInco2'
infile = '/disk/kea/WAIS/orig/xlob/ICP3/JKB2d/F56T01c/RADnh3/bxds'
unfoc(infile=infile, outdir='/tmp/outdir', stackdepth=10, incodepth=5,
      channels=channels, processes=2)

```

## Using readers for data ingest

The `unfoc` module provides several ways to read raw radar data from
the `RADnh3` and `RADnh5` data streams.  

See [docs/sphinx/source/usage_docs][usage_info] for more information.


The RadBxds class automatically detects `RADnh3` and `RADnh5` streams.  
A separate class `RADjh1Bxds` can be used with RADjh1 data.  At some point
these classes will be combined to seamlessly detect data type.

The `RadBxds` and `RADjh1Bxds` classes are the simplest to use, but
they require some preprocessing to read the radargram.


# Expected Input Data

Expected input data is a bxds file of `RADnh3`, `RADnh5`, or `bxds1/bxds2` `RADjh1` radar system
that has come from [breakout_elsa](https://github.com/UTIG/breakout_elsa) or similar.

Each of these bxds files must have an associated `ct` text file or `ct.gz` gzip-compressed
text file with ELSA `ct` metadata (timestamps and sequence numbers)


# Output Data Format

Output data files after pik1 pulse compression consist of binary array data with 4-byte
big-endian signed integer samples.  Output files include:

- `MagLoResIncoN` - binary array of radar magnitude data, dB scaled. By convention 1000*power in dB (if default `mag_scale_factor=20_000` is used)
- `PhsLoResIncoN` - binary array of radar phase data, in units of radians * `2^24`.
- `TraceNumbersN` - Text list of first ELSA sequence number used in each pulse-compressed radar trace in the binary arrays.

MagLoResIncoN and PhsLoResIncoN typically have 3200 samples per trace.

Semantic meaning of channel numbers are discussed in `doc/Multipol Data Format.pdf`.

# Acknowledgments

Supported by **NSF award 2127606** (“Collaborative Research: EarthCube Capabilities: Open Polar Radar (OPoRa) Software and Service”)
and the **G. Unger Vetlesen Antarctic Aerogeophysical Data Analysis Project** ([G. Unger Vetlesen Foundation](https://www.vetlesenfoundation.org/)).
