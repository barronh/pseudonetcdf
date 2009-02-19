__doc__ = """
.. _vertical_diffusivity
:mod:`vertical_diffusivity` -- Vertical Diffusivity File Interfaces
===================================================================

.. module:: vertical_diffusivity
   :platform: Unix, Windows
   :synopsis: Provides :ref:`PseudoNetCDF` memory map and random access read 
   based file interfaces for CAMx vertical diffusivity files.
.. moduleauthor:: Barron Henderson <barronh@unc.edu>
"""
__all__=['Memmap','Read','Transforms','Write']

import Memmap
import Read
import Write
import Transforms

if __name__ == '__main__':
    from PseudoNetCDF.camxfiles.vertical_diffusivity.Memmap import vertical_diffusivity
    from PseudoNetCDF.pncdump import pncdump_parser, \
                                    dump_from_cmd_line
    parser = pncdump_parser()
    parser.add_argument("cols", int)
    parser.add_argument("rows", int)
    (file_path, options, extra_args_dict) = parser.parse_args()

    dump_from_cmd_line(file_path, options, lambda path: vertical_diffusivity(path, **extra_args_dict))
