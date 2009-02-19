__doc__ = """
.. _wind
:mod:`wind` -- Wind File Interfaces
===================================

.. module:: wind
   :platform: Unix, Windows
   :synopsis: Provides :ref:`PseudoNetCDF` memory map and random access read 
   based file interfaces for CAMx wind files.
.. moduleauthor:: Barron Henderson <barronh@unc.edu>
"""
__all__=['Memmap','Read','Write','Transforms']

import Memmap
import Read
import Write
import Transforms

if __name__ == '__main__':
    from PseudoNetCDF.camxfiles.wind.Memmap import wind
    from PseudoNetCDF.pncdump import pncdump_parser, \
                                    dump_from_cmd_line
    parser = pncdump_parser()
    parser.add_argument("cols", int)
    parser.add_argument("rows", int)
    (file_path, options, extra_args_dict) = parser.parse_args()

    dump_from_cmd_line(file_path, options, lambda path: wind(path, **extra_args_dict))
