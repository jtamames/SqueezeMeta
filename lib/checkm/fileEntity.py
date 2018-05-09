#!/usr/bin/env python
###############################################################################
#                                                                             #
#    fileEntity.py                                                            #
#                                                                             #
#    Represent a file / folder with path etc                                  #
#                                                                             #
#    Copyright (C) Michael Imelfort                                           #
#                                                                             #
###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

__author__ = "Michael Imelfort"
__copyright__ = "Copyright 2014"
__credits__ = ["Michael Imelfort"]
__license__ = "GPLv3"
__version__ = "0.1.0"
__maintainer__ = "Michael Imelfort"
__email__ = "mike@mikeimelfort.com"
__status__ = "Beta"

###############################################################################
###############################################################################
###############################################################################
###############################################################################

# system includes
import sys
import os

# local includes

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class FileEntity(object):
    """Basic file entity"""
    def __init__(self,
                 name,      # the name of the entity on the file system ( Full path to root dir if id: ROOT)
                 path,      # the local path to this entity
                 parent,    # the entity (type == 'dir') that contains this. ( None for id: ROOT )
                 hashd,     # hash of the entity if type == 'file'
                 size       # size of the file in bytes
                 ):
        self.name = name
        self.path=path
        self.parent = parent
        self.hashd = hashd
        self.size = size

    def getFullPath(self):
        """get the full path to this entity"""
        if self.parent == None:
            return ""
        else:
            return os.path.join(self.parent.getFullPath(), self.name)

    def checkIntegrity(self):
        """Check the file for corruption"""
        if self.type == 'dir':
            return True
        else:
            # check the hashd and compare against the recorded MD5
            return True

    def __str__(self):
        if self.parent is not None:
            return "\t".join([os.path.join(self.path,self.name),self.hashd,str(self.size)])
        return ""

#------------------------------------------------------------------------------
# Handling IDs




###############################################################################
###############################################################################
###############################################################################
###############################################################################

