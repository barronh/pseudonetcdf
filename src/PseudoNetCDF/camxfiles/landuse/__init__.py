__doc__ = """
.. _landuse
:mod:`landuse` -- Landuse File Interfaces
=========================================

.. module:: landuse
   :platform: Unix, Windows
   :synopsis: Provides :ref:`PseudoNetCDF` memory map and random access read 
   based file interfaces for CAMx landuse files.
.. moduleauthor:: Barron Henderson <barronh@unc.edu>
"""
__all__=['Memmap']

import Memmap

if __name__ == '__main__':
    from PseudoNetCDF.camxfiles.landuse.Memmap import landuse
    from PseudoNetCDF.pncdump import pncdump_parser, \
                                    dump_from_cmd_line
    parser = pncdump_parser()
    parser.add_argument("cols", int)
    parser.add_argument("rows", int)
    (file_path, options, extra_args_dict) = parser.parse_args()

    dump_from_cmd_line(file_path, options, lambda path: landuse(path, **extra_args_dict))
