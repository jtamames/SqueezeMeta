import pkg_resources  # part of setuptools
__version__ = pkg_resources.require("concoct")[0].version
