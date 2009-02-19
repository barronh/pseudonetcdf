__doc__ = """
.. _humidity
:mod:`humidity` -- Humidity File Interfaces
===========================================

.. module:: humidity
   :platform: Unix, Windows
   :synopsis: Provides :ref:`PseudoNetCDF` memory map and random access read 
   based file interfaces for CAMx humidity files.
.. moduleauthor:: Barron Henderson <barronh@unc.edu>
"""
__all__=['Memmap','Read','Transforms','Write']

import Memmap
import Read
import Write
import Transforms

if __name__ == '__main__':
    from PseudoNetCDF.camxfiles.humidity.Memmap import humidity
    from PseudoNetCDF.pncdump import pncdump_parser, \
                                    dump_from_cmd_line
    parser = pncdump_parser()
    parser.add_argument("cols", int)
    parser.add_argument("rows", int)
    (file_path, options, extra_args_dict) = parser.parse_args()

    dump_from_cmd_line(file_path, options, humidity)
