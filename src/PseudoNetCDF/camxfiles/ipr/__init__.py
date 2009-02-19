__doc__ = """
.. _ipr
:mod:`ipr` -- IPR File Interfaces
=================================

.. module:: ipr
   :platform: Unix, Windows
   :synopsis: Provides :ref:`PseudoNetCDF` memory map and random access read 
   based file interfaces for CAMx ipr files.
.. moduleauthor:: Barron Henderson <barronh@unc.edu>
"""
__all__=['Memmap','Read']

import Memmap
import Read

if __name__ == '__main__':
    from PseudoNetCDF.camxfiles.ipr.Memmap import ipr
    from PseudoNetCDF.pncdump import pncdump_parser, \
                                    dump_from_cmd_line
    parser = pncdump_parser()
    (file_path, options, extra_args_dict) = parser.parse_args()

    dump_from_cmd_line(file_path, options, ipr)
