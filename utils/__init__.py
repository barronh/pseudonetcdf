__doc__ = """
.. _utils
:mod:`utils` -- pyPA utilities
==============================

.. module:: utils
   :platform: Unix, Windows
   :synopsis: Provides general utilities used in :ref:`pyPA`
.. moduleauthor:: Barron Henderson <barronh@unc.edu>
"""

__all__=['CAMxFiles', \
         'CAMxTransforms', \
         'util', \
         'xml', \
         'spatial']
         
HeadURL="$HeadURL$"
ChangeDate = "$LastChangedDate$"
RevisionNum= "$LastChangedRevision$"
ChangedBy  = "$LastChangedBy$"
__version__ = RevisionNum

import CAMxFiles
import CAMxTransforms
import util
import xml
import spatial