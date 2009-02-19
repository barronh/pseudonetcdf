__doc__ = """
.. _finst
:mod:`finst` -- Fine UAM-IV File Interfaces
===========================================

.. module:: finst
   :platform: Unix, Windows
   :synopsis: Provides :ref:`PseudoNetCDF` memory map and random access read 
   based file interfaces for CAMx fine UAM-IV files.
.. moduleauthor:: Barron Henderson <barronh@unc.edu>
"""
__all__=['Memmap']

import Memmap


if __name__ == '__main__':
    from PseudoNetCDF.camxfiles.finst.Memmap import finst
    from PseudoNetCDF.pncdump import pncdump_parser, \
                                    dump_from_cmd_line
    parser = pncdump_parser()
    parser.add_argument("cols", int)
    parser.add_argument("rows", int)
    (file_path, options, extra_args_dict) = parser.parse_args()

    dump_from_cmd_line(file_path, options, lambda path: finst)
