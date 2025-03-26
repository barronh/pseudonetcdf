__all__ = ['aermod_plotfile']

from ._aermod_plotfile import reader
from ..register import registerreader
aermod_plotfile = reader
registerreader('aermod', reader)
