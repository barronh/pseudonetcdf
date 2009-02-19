__doc__ = """
.. _temperature
:mod:`temperature` -- Temperature File Interfaces
=================================================

.. module:: temperature
   :platform: Unix, Windows
   :synopsis: Provides :ref:`PseudoNetCDF` memory map and random access read 
   based file interfaces for CAMx temperature files.
.. moduleauthor:: Barron Henderson <barronh@unc.edu>
"""
__all__=['Memmap','Read','Transforms']

import Memmap
import Read
import Transforms

if __name__ == '__main__':
    from PseudoNetCDF.camxfiles.temperature.Memmap import temperature
    from PseudoNetCDF.pncdump import pncdump_parser, \
                                    dump_from_cmd_line
    parser = pncdump_parser()
    parser.add_argument("cols", int)
    parser.add_argument("rows", int)
    (file_path, options, extra_args_dict) = parser.parse_args()

    dump_from_cmd_line(file_path, options, lambda path: temperature(path, **extra_args_dict))
