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

import time
import logging
from functools import reduce


class TimeKeeper:
    """Helper class for tracking and reporting elapsed time."""

    def __init__(self):
        """Initialization."""
        self.logger = logging.getLogger()
        self.start_time = time.time()
        self.last_log_time = self.start_time

    def start_timer(self):
        """Restart the timer"""
        self.start_time = time.time()
        self.last_log_time = self.start_time

    def get_time_stamp(self):
        """Make a time stamp"""
        now = time.time()
        ret_str = "\n  { Current stage: %s || Total: %s }" % (self.seconds_to_str(now - self.last_log_time), self.seconds_to_str(now - self.start_time))
        self.last_log_time = now
        return ret_str

    def print_time_stamp(self):
        """Print time stamp."""
        self.logger.info(self.get_time_stamp())

    def seconds_to_str(self, time_in_seconds):
        """Convert elapsed time in seconds to human readable string.

        Parameters
        ----------
        time_in_seconds
            Time in seconds.
        """
        rediv = lambda ll, b: list(divmod(ll[0], b)) + ll[1:]
        return "%d:%02d:%02d.%03d" % tuple(reduce(rediv, [[time_in_seconds * 1000, ], 1000, 60, 60]))
