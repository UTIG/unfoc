.. unfoc documentation master file

unfoc
=====

Welcome to **unfoc** — an unfocused radar processing package for ice-penetrating
radar data collected by the UTIG/SOAR systems.

The primary goal is to load raw radar traces into NumPy arrays together with
their timing/sequence metadata (CT) and then produce a pulse-compressed
unfocused product (often called *pik1*).

Features
--------

- Reading radar data
  - RADnh3, RADnh5, RADjh1 
- Metadata association via CT files
- Noise filtering and denoising (burst suppression)
- Channel parsing and trace stacking
- Writing unfocused radargrams (magnitude/phase/trace numbers)

Project links
-------------

- GitHub repository: https://github.com/UTIG/unfoc
- Issue tracker: https://github.com/UTIG/unfoc/issues

Quick start
-----------

Install and run (see detailed steps in :doc:`installation` and :doc:`usage_guide`):

.. code-block:: bash

   # clone (private org access required) or install once published to PyPI
   git clone https://github.com/UTIG/unfoc.git
   cd unfoc
   pip install .

   # process a bxds file to an output directory
   run_unfoc.py -i /path/to/bxds -o /tmp/outdir \
                --stackdepth 10 --incodepth 5 \
                --channels LoResInco1,LoResInco2


References
----------

- Peters, M. E., Blankenship, D. D., & Morse, D. L. (2005).
  *Analysis techniques for coherent airborne radar sounding: Application to West Antarctic ice streams*.
  **Journal of Geophysical Research: Solid Earth**, 110(B6), B06303.
  https://doi.org/10.1029/2004JB003222

Funding
-------

This work has been supported by:


- NSF Award 2127606 — “Collaborative Research: EarthCube Capabilities: Open Polar Radar (OPoRa) Software and Service”.
- G. Unger Vetlesen Antarctic Aerogeophysical Data Analysis Project (`G. Unger Vetlesen Foundation <https://www.vetlesenfoundation.org/>`_).

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   usage_guide
   package_design
   modules

Indices and tables
------------------

* :ref:`modindex`
