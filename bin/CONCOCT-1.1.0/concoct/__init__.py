#import pkg_resources  # part of setuptools
#__version__ = pkg_resources.require("concoct")[0].version
from importlib import metadata as importlib_metadata
__version__ = importlib_metadata.version('concoct')
