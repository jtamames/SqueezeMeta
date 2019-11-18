import os

# get package version
packageDir = os.path.dirname(os.path.realpath(__file__))
versionFile = open(os.path.join(packageDir, 'VERSION'))
__version__ = versionFile.read().strip()
