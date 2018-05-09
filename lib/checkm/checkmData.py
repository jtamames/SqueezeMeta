###############################################################################
#
# checkmData.py - database management utilities
#
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

import os
import sys
import logging
from pkg_resources import resource_filename
import json

import checkm.manifestManager as mm


class DBConfig(object):
    """CheckM uses packageutils to distribute a file called "DATA_CONFIG" which is placed with
    the other files during installation. This file stores information about where data is
    stored locally and where to look for remote updates. This class essentially exposes
    the DATA_CONFIG file as an object.
    """
    def __init__(self):
        self.logger = logging.getLogger()
        self.configFile = os.path.abspath(resource_filename('checkm', 'DATA_CONFIG'))
        self.values = self.getConfig()

#-----------------------------------------------------------------------------
# Read and write local config file

    def getConfig(self):
        """Get a listing of the versions of files in the local data config"""
        try:
            with open(self.configFile, 'r') as local_config:
                # config is a one line file
                for line in local_config:
                    return json.loads(line)

        except Exception:
            self.logger.error("There seems to be a problem with loading the CheckM config file")
            self.logger.error("Please check the permissions / existence / contents of:")
            self.logger.error(self.configFile)
            raise

        return {}

    def setConfig(self):
        """Update the local config to reflect and changes made"""
        if self.checkPermissions():
            with open(self.configFile, 'w') as config_fh:
                config_fh.write(json.dumps(self.values))

#-----------------------------------------------------------------------------
# Display config

    def displayConfig(self):
        """Print out the contents of the CheckM config file"""
        self.logger.info("Contents of the config located at:\n\n\t%s\n" % (self.configFile))
        self.logger.info("-------------------------------------------------------------------------------")
        self.logger.info("Data root: %s" % self.values["dataRoot"])
        self.logger.info("Local manifest name: %s" % self.values["remoteManifestURL"])
        self.logger.info("Manifest type: %s" % self.values["remoteManifestName"])
        self.logger.info("Remote URL: %s" % self.values["localManifestName"])
        self.logger.info("Remote manifest name: %s" % self.values["manifestType"])
        self.logger.info("-------------------------------------------------------------------------------")

#-----------------------------------------------------------------------------
# Housekeeping

    def checkPermissions(self):
        """Work out if we have permission to write to the CheckM config before attempting to make changes"""
        try:
            open(self.configFile, 'a')
        except IOError, e:
            print "You do not seem to have permission to edit the checkm config file"
            print "located at %s" % self.configFile
            print "Please try again with updated privileges. Error was:\n"
            print e
            return False
        return True


class DBManager(mm.ManifestManager):

    """Manage all aspects of data location and version control."""
    def __init__(self):
        mm.ManifestManager.__init__(self, timeout=15)
        self.logger = logging.getLogger()
        self.config = DBConfig()  # load inbuilt configuration
        self.type = self.config.values["manifestType"]

        # check that the data root is legit
        manifestFile = os.path.join(self.config.values["dataRoot"], mm.__MANIFEST__)
        if not os.path.exists(self.config.values["dataRoot"]) or not os.path.exists(manifestFile):
            self.config.values["dataRoot"] = ""

        if self.config.values["dataRoot"] == "":
            # no data folder set.
            print "It seems that the CheckM data folder has not been set yet or has been removed. Running: 'checkm data setRoot'."
            if not self.setRoot():
                print "Sorry, CheckM cannot run without a valid data folder."

    def runAction(self, action):
        """Main entry point for the updating code"""
        
        if action[0] == "setRoot":
            if len(action) > 1:
                path = self.setRoot(path=action[1])
            else:
                path = self.setRoot()

            if path is None:
                self.logger.info("Data location not changed")
            else:
                self.logger.info("Data location successfully changed to: %s" % path)
        else:
            self.logger.error("Unknown action: %s" % action[0])

    def setRoot(self, path=None):
        """Set the data folder"""
        # check to see we have permission to update the data root
        if not self.config.checkPermissions():
            return None

        # validate the supplied path or prompt for a new one
        path = self.confirmPath(path=path)
        if path is None:
            # The user is not interested
            return None

        # path should be set, exist and be writable
        self.config.values["dataRoot"] = path

        # save the new path
        self.config.setConfig()

        return path

    def confirmPath(self, path=None):
        """Ask the user to supply a path"""
        path_set = False
        minimal = False
        while not path_set:
            if not path:
                if(minimal):
                    path = raw_input("Please specify a location or type 'abort' to stop trying: \n")
                else:
                    path = raw_input("Where should CheckM store it's data?\n" \
                                    "Please specify a location or type 'abort' to stop trying: \n")

            if path.upper() == "ABORT":
                # user has given up
                return None
            else:
                path = os.path.abspath(os.path.expanduser(path))

            print ""
            if os.path.exists(path):
                # path exists
                if os.access(path, os.W_OK):
                    # path is writable
                    path_set = True
                    print "Path [%s] exists and you have permission to write to this folder." % path
                else:
                    print "Path [%s] exists but you do not have permission to write to this folder." % path
            else:
                # path does not exist, try to make it
                "Path [%s] does not exist so I will attempt to create it" % path
                try:
                    self.makeSurePathExists(path)
                    print "Path [%s] has been created and you have permission to write to this folder." % path
                    path_set = True
                except Exception:
                    print "Unable to make the folder, Error was: %s" % sys.exc_info()[0]
                minimal = True

        # (re)make the manifest file
        print "(re) creating manifest file (please be patient)."
        self.createManifest(path, self.config.values["localManifestName"])

        return path

    def checkPermissions(self):
        """See if the user has permission to write to the data directory"""
        if not os.access(self.config.values["dataRoot"], os.W_OK):
            print "You do not seem to have permission to edit the CheckM data folder"
            print "located at %s" % self.config.values["dataRoot"]
            return False

        return True
