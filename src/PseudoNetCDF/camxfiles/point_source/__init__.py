__doc__ = """
.. _point_source
:mod:`point_source` -- Point Source File Interfaces
===================================================

.. module:: point_source
   :platform: Unix, Windows
   :synopsis: Provides :ref:`PseudoNetCDF` memory map and random access read 
   based file interfaces for CAMx point source files.
.. moduleauthor:: Barron Henderson <barronh@unc.edu>
"""
__all__=['Memmap','Read','Write']

import Memmap
import Read
import Write

if __name__ == '__main__':
    from PseudoNetCDF.camxfiles.point_source.Memmap import point_source
    from PseudoNetCDF.pncdump import pncdump_parser, \
                                    dump_from_cmd_line
    parser = pncdump_parser()
    (file_path, options, extra_args_dict) = parser.parse_args()

    dump_from_cmd_line(file_path, options, point_source)
