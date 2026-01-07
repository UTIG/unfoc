Installation
============

`utig-unfoc` is the Unfocused Radar Processor developed at the University of Texas Institute for Geophysics (UTIG).
It is currently in **beta** and available in the UTIG GitHub organization:
https://github.com/UTIG/unfoc

You must be added to the **UTIG GitHub organization** to access the repository. Please contact the maintainers for access.

Once added, you can install the package in one of the following ways:

--------------------------------
1. Manual Installation (Release)
--------------------------------

To manually install from the official release (e.g., v2.1.0):

**Option A – Install from `.whl` file:**

1. Download the wheel file (`utig_unfoc-2.1.0-py3-none-any.whl`) from:
   https://github.com/UTIG/unfoc/releases/tag/v2.1.0

2. Run:

.. code-block:: bash

    pip install utig_unfoc-2.1.0-py3-none-any.whl

**Option B – Install from source tarball:**

.. code-block:: bash

    tar -xzf utig_unfoc-2.1.0.tar.gz
    cd utig_unfoc-2.1.0
    pip install .

--------------------------------
2. Clone and Install from GitHub
--------------------------------

If you have access to the private GitHub repository (you must be added to the
``UTIG`` organization first), you can install the **latest** version in a few ways.

.. code-block:: bash

    git clone git@github.com:UTIG/unfoc.git
    cd unfoc
    pip install .

You may also use HTTPS:

.. code-block:: bash

    git clone https://github.com/UTIG/unfoc.git

----------------------------------
3. PyPI Installation (Coming Soon)
----------------------------------

Once published to PyPI, you will be able to install with:

.. code-block:: bash

    pip install utexas-ipr-unfoc

Stay tuned for announcements on PyPI availability.

------------
Requirements
------------

**Python version:** 3.7–3.12  
**Dependencies:**

- numpy ≥ 1.19.4
- scipy ≥ 1.5.4

Optional development dependencies:

- coverage ≥ 7.2.7 (for test coverage reporting)


