__all__ = ['version']

try:
    from importlib.metadata import version as _version
except ImportError:
    # if the fallback library is missing, version will be incorrectly
    # reported as 0.0.0
    from importlib_metadata import version as _version

try:
    version = _version("PseudoNetCDF")
except Exception:
    # Local copy or not installed with setuptools.
    # Disable minimum version checks on downstream libraries.
    version = "9.9.9"
