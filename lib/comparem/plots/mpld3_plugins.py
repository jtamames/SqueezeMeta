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

import matplotlib
import mpld3


class LinkedBrush(mpld3.plugins.PluginBase):
    JAVASCRIPT = """
        mpld3.LinkedBrushPlugin = refinem_LinkedBrushPlugin;
        mpld3.register_plugin("refinem_linkedbrush", refinem_LinkedBrushPlugin);
        refinem_LinkedBrushPlugin.prototype = Object.create(mpld3.Plugin.prototype);
        refinem_LinkedBrushPlugin.prototype.constructor = refinem_LinkedBrushPlugin;
        refinem_LinkedBrushPlugin.prototype.requiredProps = [ "id" ];
        refinem_LinkedBrushPlugin.prototype.defaultProps = {
          button: true,
          enabled: null
        };

        function refinem_LinkedBrushPlugin(fig, props) {
          mpld3.Plugin.call(this, fig, props);
          if (this.props.enabled === null) {
            this.props.enabled = !this.props.button;
          }
          var enabled = this.props.enabled;
          if (this.props.button) {
            var BrushButton = mpld3.ButtonFactory({
              buttonID: "refinem_linkedbrush",
              sticky: true,
              actions: [ "drag" ],
              onActivate: this.activate.bind(this),
              onDeactivate: this.deactivate.bind(this),
              onDraw: function() {
                this.setState(enabled);
              },
              icon: function() {
                return mpld3.icons["brush"];
              }
            });
            this.fig.buttons.push(BrushButton);
          }
          this.extentClass = "refinem_linkedbrush";
        }

        refinem_LinkedBrushPlugin.prototype.activate = function() {
          if (this.enable) this.enable();
        };

        refinem_LinkedBrushPlugin.prototype.deactivate = function() {
          if (this.disable) this.disable();
        };

        refinem_LinkedBrushPlugin.prototype.get_selected = function() {
          if (this.get_selected) this.get_selected();
        };

        refinem_LinkedBrushPlugin.prototype.draw = function() {
          var obj = mpld3.get_element(this.props.id);
          if (obj === null) {
            throw "LinkedBrush: no object with id='" + this.props.id + "' was found";
          }

          var fig = this.fig;
          if (!("offsets" in obj.props)) {
            throw "Plot object with id='" + this.props.id + "' is not a scatter plot";
          }

          var dataKey = "offsets" in obj.props ? "offsets" : "data";
          mpld3.insert_css("#" + fig.figid + " rect.extent." + this.extentClass, {
            fill: "#000",
            "fill-opacity": .05,
            stroke: "#fff"
          });

          mpld3.insert_css("#" + fig.figid + " path.mpld3-unselected", {
            opacity: .2
          });

          var dataClass = "mpld3data-" + obj.props[dataKey];
          var brush = fig.getBrush();
          var dataByAx = [];
          fig.axes.forEach(function(ax) {
            var axData = [];
            ax.elements.forEach(function(el) {
              if (el.props[dataKey] === obj.props[dataKey]) {
                el.group.classed(dataClass, true);
                axData.push(el);
              }
            });
            dataByAx.push(axData);
          });

          var allData = [];
          var selectedData = fig.canvas.selectAll("." + dataClass);
          var unselectedData = fig.canvas.selectAll("." + dataClass);
          var currentAxes;
          function brushstart(d) {
            if (currentAxes != this) {
              d3.select(currentAxes).call(brush.clear());
              currentAxes = this;
              brush.x(d.xdom).y(d.ydom);
            }
          }

          function brushmove(d) {
            var data = dataByAx[d.axnum];
            if (data.length > 0) {
              var ix = data[0].props.xindex;
              var iy = data[0].props.yindex;
              var e = brush.extent();

              if (brush.empty()) {
                selectedData.selectAll("path").classed("mpld3-selected", false);
                unselectedData.selectAll("path").classed("mpld3-unselected", false);
              } else {
                selectedData.selectAll("path").classed("mpld3-selected", function(p) {
                  return !(e[0][0] > p[ix] || e[1][0] < p[ix] || e[0][1] > p[iy] || e[1][1] < p[iy]);
                });

                unselectedData.selectAll("path").classed("mpld3-unselected", function(p) {
                  return e[0][0] > p[ix] || e[1][0] < p[ix] || e[0][1] > p[iy] || e[1][1] < p[iy];
                });
              }
            }
          }

          function brushend(d) {
            if (brush.empty()) {
              selectedData.selectAll("path").classed("mpld3-selected", false);
              unselectedData.selectAll("path").classed("mpld3-selected", false);
            }
          }

          this.enable = function() {
            this.fig.showBrush(this.extentClass);
            brush.on("brushstart", brushstart).on("brush", brushmove).on("brushend", brushend);
            this.enabled = true;
          };

          this.disable = function() {
            d3.select(currentAxes).call(brush.clear());
            this.fig.hideBrush(this.extentClass);
            this.enabled = false;
          };

          this.disable();
        };
    """

    def __init__(self, points, button=True, enabled=True):
        if isinstance(points, matplotlib.lines.Line2D):
            suffix = "pts"
        else:
            suffix = None
        self.dict_ = {"type": "refinem_linkedbrush",
                      "button": button,
                      "enabled": False,
                      "id": mpld3.utils.get_id(points, suffix)}


class Tooltip(mpld3.plugins.PluginBase):
    """A Plugin to enable an HTML tooltip.

    This extends the PointHTMLTooltip class in mpld3. It adds
    a mousedown() event which writes the label of clicked
    points to an HTML element with the id 'selected_points'.

    formated text which hovers over points.
    Parameters
    ----------
    points : matplotlib Collection or Line2D object
        The figure element to apply the tooltip to
    labels : list
        The labels for each point in points, as strings of unescaped HTML.
    hoffset, voffset : integer, optional
        The number of pixels to offset the tooltip text.  Default is
        hoffset = 0, voffset = 10
    css : str, optional
        css to be included, for styling the label html if desired
    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> from mpld3 import fig_to_html, plugins
    >>> fig, ax = plt.subplots()
    >>> points = ax.plot(range(10), 'o')
    >>> labels = ['<h1>{title}</h1>'.format(title=i) for i in range(10)]
    >>> plugins.connect(fig, Tooltip(points[0], labels))
    >>> fig_to_html(fig)
    """

    JAVASCRIPT = """
        mpld3.register_plugin("tooltip", Tooltip);
        Tooltip.prototype = Object.create(mpld3.Plugin.prototype);
        Tooltip.prototype.constructor = Tooltip;
        Tooltip.prototype.requiredProps = ["id"];
        Tooltip.prototype.defaultProps = {labels:null,
                                                    hoffset:0,
                                                    voffset:10};
        function Tooltip(fig, props){
            mpld3.Plugin.call(this, fig, props);
        };
        Tooltip.prototype.draw = function(){
           var obj = mpld3.get_element(this.props.id);
           var labels = this.props.labels;
           var tooltip = d3.select("body").append("div")
                        .attr("class", "mpld3-tooltip")
                        .style("position", "absolute")
                        .style("z-index", "10")
                        .style("visibility", "hidden");
           obj.elements()
               .on("mousedown", function(d, i){
                                var div = document.getElementById("selected_points");
                                selected_points[selected_points.length] = labels[i]
                                div.innerHTML = selected_points.join('<br>');})
               .on("mouseover", function(d, i){
                                  tooltip.html(labels[i])
                                         .style("visibility", "visible");})
               .on("mousemove", function(d, i){
                      tooltip
                        .style("top", d3.event.pageY + this.props.voffset + "px")
                        .style("left",d3.event.pageX + this.props.hoffset + "px");
                     }.bind(this))
               .on("mouseout",  function(d, i){
                               tooltip.style("visibility", "hidden");});
        };
    """

    def __init__(self, points, labels=None,
                 hoffset=0, voffset=10, css=None):
        self.points = points
        self.labels = labels
        self.voffset = voffset
        self.hoffset = hoffset
        self.css_ = css or ""
        if isinstance(points, matplotlib.lines.Line2D):
            suffix = "pts"
        else:
            suffix = None
        self.dict_ = {"type": "tooltip",
                      "id": mpld3.utils.get_id(points, suffix),
                      "labels": labels,
                      "hoffset": hoffset,
                      "voffset": voffset}

    # Additional script to add at global scope.
    script_global = ('var selected_points = []\n'
                        'function clear_selection_list() {\n'
                            'selected_points = []\n'
                            'var div = document.getElementById("selected_points");\n'
                            'div.innerHTML = "";\n'
                        '};\n')

    # Additional HTML to add to body
    html_body = ('<hr>\n'
                    '<div id="selected_points"></div>\n'
                    '<br>\n'
                    '<button onclick="clear_selection_list()">Clear</button>\n')
