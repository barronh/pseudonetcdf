__doc__ = """
.. _height_pressure
:mod:`height_pressure` -- Height/Pressure File Interfaces
=========================================================

.. module:: height_pressure
   :platform: Unix, Windows
   :synopsis: Provides :ref:`PseudoNetCDF` memory map and random access read 
   based file interfaces for CAMx height pressure files.
.. moduleauthor:: Barron Henderson <barronh@unc.edu>
"""
__all__=['Memmap','Read','Write','Transforms']

import Memmap
import Read
import Write
import Transforms

if __name__ == '__main__':
    from PseudoNetCDF.camxfiles.height_pressure.Memmap import height_pressure
    from PseudoNetCDF.pncdump import pncdump_parser, \
                                    dump_from_cmd_line
    parser = pncdump_parser()
    parser.add_argument("cols", int)
    parser.add_argument("rows", int)
    (file_path, options, extra_args_dict) = parser.parse_args()

    dump_from_cmd_line(file_path, options, lambda path: height_pressure(path, **extra_args_dict))
