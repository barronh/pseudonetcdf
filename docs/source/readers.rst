.. Readers

Formats
~~~~~~~

PseudoNetCDF helps make files act like netCDF by making readers that convert the
binary files data data and metadata, and by providing convenience functions for
the metadata. To do that, it has many `format` options. The easiest way to use
them is via :func:`pncopen`.

.. code-block:: python

    import PseudoNetCDF as pnc
    path = '/path/to/file'
    fmt = 'ioapi'  # ioapi is for CMAQ; bpch is for GEOS-Chem; CAMx has options...
    f = pnc.pncopen(path, format='?')


You can find a full list of options using `pncopen(help=True)`. The list below
gives common `format` options for specific types of files.

.. toctree::
   :maxdepth: 1
   :glob:

   readers/*

