"""
.. _cloud_rain
:mod:`cloud_rain` -- Cloud/Rain File Interfaces
===============================================

.. module:: cloud_rain
   :platform: Unix, Windows
   :synopsis: Provides :ref:`PseudoNetCDF` file interfaces for CAMx cloud/rain files.
.. moduleauthor:: Barron Henderson <barronh@unc.edu>
"""
__all__=['Memmap','Transforms']

import Memmap
import Transforms

if __name__ == '__main__':
    from PseudoNetCDF.camxfiles.cloud_rain.Memmap import cloud_rain
    from PseudoNetCDF.pncdump import pncdump_parser, \
                                    dump_from_cmd_line
    parser = pncdump_parser()
    (file_path, options, extra_args_dict) = parser.parse_args()

    dump_from_cmd_line(file_path, options, cloud_rain)
