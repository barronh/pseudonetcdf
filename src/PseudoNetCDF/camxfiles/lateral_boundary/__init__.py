__doc__ = """
.. _uamiv
:mod:`uamiv` -- UAM-IV File Interfaces
======================================

.. module:: uamiv
   :platform: Unix, Windows
   :synopsis: Provides :ref:`PseudoNetCDF` memory map and random access read 
   based file interfaces for CAMx UAM-IV files.
.. moduleauthor:: Barron Henderson <barronh@unc.edu>
"""
__all__=['Memmap','Read','Write','Transforms']

import Memmap

# _camx_units is based on file name and an aerosol flag (True = aerosol)
if __name__ == '__main__':
    from PseudoNetCDF.camxfiles.lateral_boundary.Memmap import lateral_boundary
    from PseudoNetCDF.pncdump import pncdump_parser, \
                                    dump_from_cmd_line
    parser = pncdump_parser()
    (file_path, options, extra_args_dict) = parser.parse_args()

    dump_from_cmd_line(file_path, options, lateral_boundary)
