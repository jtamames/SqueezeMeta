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

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2014'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'

import sys
from collections import namedtuple

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.transforms as mtransforms

from matplotlib.patches import Rectangle

import matplotlib as mpl

import numpy as np


class AbstractPlot(FigureCanvas):
    """Abstract base class for plotting."""
    
    Options = namedtuple('Options', 'width height label_font_size tick_font_size dpi')

    def __init__(self, options):
        """Initialization."""

        if options:
            self.options = options
        else:
            Options = namedtuple('Options', 'width height label_font_size tick_font_size dpi')
            self.options = Options(6, 6, 10, 8, 300)

        self.set_font_size(self.options.label_font_size, self.options.tick_font_size)
        self.fig = Figure(facecolor='white', dpi=self.options.dpi)

        FigureCanvas.__init__(self, self.fig)

        self.cid = None

        self.type = '<none>'
        self.name = '<none>'

        self.axes_colour = (0.5, 0.5, 0.5)

    def set_font_size(self, label_size, tick_size):
        """Set font size for all text elements."""
        mpl.rcParams['font.size'] = label_size
        mpl.rcParams['axes.titlesize'] = label_size
        mpl.rcParams['axes.labelsize'] = label_size
        mpl.rcParams['xtick.labelsize'] = tick_size
        mpl.rcParams['ytick.labelsize'] = tick_size
        mpl.rcParams['legend.fontsize'] = label_size
        mpl.rcParams['svg.fonttype'] = 'none'

    def save_plot(self, filename, dpi=300):
        """Save plot to file."""

        imgFormat = filename[filename.rfind('.') + 1:len(filename)]
        if imgFormat in ['png', 'pdf', 'ps', 'eps', 'svg']:
            self.fig.savefig(filename, format=imgFormat, dpi=dpi, facecolor='white', edgecolor='white', bbox_inches='tight')
        else:
            pass

    def save_html(self, output_html, html_script=None, html_body=None):
        """Save figure as HTML.

        Parameters
        ----------
        output_html : str
            Name of output file.
        html_script : str
            Additional java script to append to script section.
        html_body : str
            Additional HTML to append to end of the body section.
        """

        try:
            __import__('mpld3')
        except ImportError:
            print('[Error] The mpld3 module is required to save HTML figures.')
            sys.exit()

        import mpld3

        # modify figure properties for better web viewing
        self.fig.dpi = 96
        html_str = mpld3.fig_to_html(self.fig, template_type='simple')

        if html_script:
            html_script_start = html_str.find('<script type="text/javascript">') + len('<script type="text/javascript">')
            html_str = html_str[0:html_script_start] + '\n' + html_script + '\n' + html_str[html_script_start:]

        if html_body:
            html_str += '\n<body>\n' + html_body + '\n</body>\n'

        html_str = '<center>' + html_str + '</center>'

        fout = open(output_html, 'w')
        fout.write(html_str)
        fout.close()

        # restore figure properties
        self.fig.dpi = self.options.dpi

    def prettify(self, axis):
        """Modify axis properties to make a cleaner plot."""

        for a in axis.yaxis.majorTicks:
            a.tick1On = True
            a.tick2On = False

        for a in axis.xaxis.majorTicks:
            a.tick1On = True
            a.tick2On = False

        for line in axis.yaxis.get_ticklines():
            line.set_color(self.axes_colour)

        for line in axis.xaxis.get_ticklines():
            line.set_color(self.axes_colour)

        for loc, spine in axis.spines.items():
            if loc in ['right', 'top']:
                spine.set_color('none')
            else:
                spine.set_color(self.axes_colour)

    def label_extents(self, xLabels, xFontSize, xRotation, yLabels, yFontSize, yRotation):
        """Get maximum extents of labels."""

        self.fig.clear()

        tempAxes = self.fig.add_axes([0, 0, 1.0, 1.0])

        tempAxes.set_xticks(np.arange(len(xLabels)))
        tempAxes.set_yticks(np.arange(len(yLabels)))

        xText = tempAxes.set_xticklabels(xLabels, size=xFontSize, rotation=xRotation)
        yText = tempAxes.set_yticklabels(yLabels, size=yFontSize, rotation=yRotation)

        bboxes = []
        for label in xText:
            bbox = label.get_window_extent(self.get_renderer())
            bboxi = bbox.inverse_transformed(self.fig.transFigure)
            bboxes.append(bboxi)
        xLabelBounds = mtransforms.Bbox.union(bboxes)

        bboxes = []
        for label in yText:
            bbox = label.get_window_extent(self.get_renderer())
            bboxi = bbox.inverse_transformed(self.fig.transFigure)
            bboxes.append(bboxi)
        yLabelBounds = mtransforms.Bbox.union(bboxes)

        self.fig.clear()

        return xLabelBounds, yLabelBounds

    def x_label_extents(self, labels, fontSize, rotation=0):
        """Get maximum extents of x-axis labels."""

        self.fig.clear()

        tempAxes = self.fig.add_axes([0, 0, 1.0, 1.0])
        tempAxes.set_xticks(np.arange(len(labels)))
        xLabels = tempAxes.set_xticklabels(labels, size=fontSize, rotation=rotation)

        bboxes = []
        for label in xLabels:
            bbox = label.get_window_extent(self.get_renderer())
            bboxi = bbox.inverse_transformed(self.fig.transFigure)
            bboxes.append(bboxi)
        xLabelBounds = mtransforms.Bbox.union(bboxes)

        self.fig.clear()

        return xLabelBounds

    def y_label_extents(self, labels, fontSize, rotation=0):
        """Get maximum extents of y-axis labels."""

        self.fig.clear()

        tempAxes = self.fig.add_axes([0, 0, 1.0, 1.0])
        tempAxes.set_yticks(np.arange(len(labels)))
        yLabels = tempAxes.set_yticklabels(labels, size=fontSize, rotation=rotation)

        bboxes = []
        for label in yLabels:
            bbox = label.get_window_extent(self.get_renderer())
            bboxi = bbox.inverse_transformed(self.fig.transFigure)
            bboxes.append(bboxi)
        yLabelBounds = mtransforms.Bbox.union(bboxes)

        self.fig.clear()

        return yLabelBounds

    def format_labels(self, labels):
        """Make labels friendly to humans."""

        formattedLabels = []
        for label in labels:
            value = float(label.get_text())
            if value < 0.01:
                valueStr = '%.2e' % value
                if 'e-00' in valueStr:
                    valueStr = valueStr.replace('e-00', 'e-')
                elif 'e-0' in valueStr:
                    valueStr = valueStr.replace('e-0', 'e-')
            else:
                valueStr = '%.3f' % value

            formattedLabels.append(valueStr)

        return formattedLabels

    def remove_extra_zeros(self, label):
        """Remove unncessary trailing zeros from labels."""

        if '.' in label:
            while label[-1] == '0':
                label = label[0:-1]

        if label[-1] == '.':  # remove potential trailing decimal point
            label = label[0:-1]

        return label

    def bounding_box(self, data, ax, label, bBoundingBoxes, bLabels):
        """Draw bounding box around data."""

        data = np.array(data)

        width = max(data[:, 0]) - min(data[:, 0])
        height = max(data[:, 1]) - min(data[:, 1])
        r = Rectangle((min(data[:, 0]), min(data[:, 1])), width, height)

        if bBoundingBoxes:
            ax.add_artist(r)
            r.set_clip_box(ax.bbox)
            r.set_alpha(0.1)
            r.set_facecolor((0.5, 0.5, 0.5))

        if bLabels:
            ax.annotate(label, xy=(min(data[:, 0]), max(data[:, 1])), xytext=(0, 0),
                            textcoords='offset points', ha='right', va='bottom',
                            bbox=dict(boxstyle='round,pad=0.5', fc=(0.5, 0.5, 0.5), alpha=0.1))
