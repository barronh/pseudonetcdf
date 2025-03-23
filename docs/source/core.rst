.. Core

Core Objects
------------

Extensive detail on the :doc:`PseudoNetCDF module<api/PseudoNetCDF>` as a whole is available.
This section documents the core that everyone uses. Even if you use CMAQ, CAMx,
GEOS-Chem, ICARTT, ... or any specific file... most of its the time you need to
know three things:

  * Question: I have a file. How do I open it?
  * Answer: Use :func:`PseudoNetCDF.pncopen` with the path and the format.

  .. code-block:: python

    import PseudoNetCDF as pnc
    path = '/path/to/file'
    fmt = 'format'  # e.g., netcdf, ioapi, ffi1001, bpch, uamiv, ...
    f = pnc.pncopen(path, format=fmt)

    # For a list of formats, run `pnc.pncopen(help=True)`


  * Question: What is `f`?
  * Answer: It is a :class:`PseudoNetCDFFile <PseudoNetCDF.core.PseudoNetCDFFile>`, which

      * has dimensions, variables, and properties.

        * The variables dictionary has :class:`PseudoNetCDFVariable  <PseudoNetCDF.core.PseudoNetCDFVariable>` objects that have:

          * dimensions, properties, and values.
          * they also act like numpy arrays.

        * The dimenions dictionary has :class:`PseudoNetCDFDimension  <PseudoNetCDF.core.PseudoNetCDFDimension>` objects that have:
          * length, and
          isunlimited property

      * The file object can be subset:

        * by variable (see :func:`subset <PseudoNetCDF.core.PseudoNetCDFFile.subsetVariables>`)
        * by dimension (applies to all variables; see :func:`slice <PseudoNetCDF.core.PseudoNetCDFFile.sliceDimensions>`)

      * The file object can apply functions to dimensions (applies to all variables; see :func:`apply <PseudoNetCDF.core.PseudoNetCDFFile.applyAlongDimensions>`)

        * mean, sum, min, max, std, ... any numpy by name, or
        * any custom function,
        * for example: `mf = f.apply(time='mean')`
      
      * The file object can be interpolated along dimensions (see :func:`interpDimension <PseudoNetCDF.core.PseudoNetCDFFile.interpDimension>`)
      * The file object can do so much more...


If you have more questions, go back to the :doc:`Example scripts <examples>` or the :doc:`API <api/PseudoNetCDF>`.
If you still have questions, open an issue on github.
