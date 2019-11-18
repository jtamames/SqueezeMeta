###############################################################################
#                                                                             #
#    color_brewer.py                                                          #
#                                                                             #
#    Easy access to the colorbrewer2.org color maps.                          #
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
__license__ = "GPL3"
__version__ = "0.2.4"
__maintainer__ = "Michael Imelfort"
__email__ = "mike@mikeimelfort.com"
__status__ = "Released"

###############################################################################


class ColorBrewer(object):
    """Provides access to the color sets specified by ColorBrewer.

        Brewer CA, Hatchard GW, and Harrower MA. 2003. ColorBrewer in Print:
          A Catalog of Color Schemes for Maps. Cartography and Geographic
          Information Science 30: 5-32.

        http://colorbrewer2.org/

        Example:

        from color_brewer import ColorBrewer
        cb = ColorBrewer()
        colors = cb.maps["qualSet1"].values()[0:3]
    """

    def __init__(self):
        """Initialize color sets."""

        self.maps = {
            "seqYellowGreen": {
                1: "#ffffe5",
                2: "#f7fcb9",
                3: "#d9f0a3",
                4: "#addd8e",
                5: "#78c679",
                6: "#41ab5d",
                7: "#238443",
                8: "#006837",
                9: "#004529"},

            "seqYellowGreenBlue": {
               1: "#ffffd9",
                2: "#edf8b1",
                3: "#c7e9b4",
                4: "#7fcdbb",
                5: "#41b6c4",
                6: "#1d91c0",
                7: "#225ea8",
                8: "#253494",
                9: "#081d58"},

            "seqGreenBlue": {
                1: "#f7fcf0",
                2: "#e0f3db",
                3: "#ccebc5",
                4: "#a8ddb5",
                5: "#7bccc4",
                6: "#4eb3d3",
                7: "#2b8cbe",
                8: "#0868ac",
                9: "#084081"},

            "seqBlueGreen": {
                1: "#f7fcfd",
                2: "#e5f5f9",
                3: "#ccece6",
                4: "#99d8c9",
                5: "#66c2a4",
                6: "#41ae76",
                7: "#238b45",
                8: "#006d2c",
                9: "#00441b"},

            "seqPurpleBlueGreen": {
                1: "#fff7fb",
                2: "#ece2f0",
                3: "#d0d1e6",
                4: "#a6bddb",
                5: "#67a9cf",
                6: "#3690c0",
                7: "#02818a",
                8: "#016c59",
                9: "#014636"},

            "seqPurpleBlue": {
                1: "#fff7fb",
                2: "#ece7f2",
                3: "#d0d1e6",
                4: "#a6bddb",
                5: "#74a9cf",
                6: "#3690c0",
                7: "#0570b0",
                8: "#045a8d",
                9: "#023858"},

            "seqBluePurple": {
                1: "#f7fcfd",
                2: "#e0ecf4",
                3: "#bfd3e6",
                4: "#9ebcda",
                5: "#8c96c6",
                6: "#8c6bb1",
                7: "#88419d",
                8: "#810f7c",
                9: "#4d004b"},

            "seqRedPurple": {
                1: "#fff7f3",
                2: "#fde0dd",
                3: "#fcc5c0",
                4: "#fa9fb5",
                5: "#f768a1",
                6: "#dd3497",
                7: "#ae017e",
                8: "#7a0177",
                9: "#49006a"},

            "seqPurpleRed": {
                1: "#f7f4f9",
                2: "#e7e1ef",
                3: "#d4b9da",
                4: "#c994c7",
                5: "#df65b0",
                6: "#e7298a",
                7: "#ce1256",
                8: "#980043",
                9: "#67001f"},

            "seqOrangeRed": {
                1: "#fff7ec",
                2: "#fee8c8",
                3: "#fdd49e",
                4: "#fdbb84",
                5: "#fc8d59",
                6: "#ef6548",
                7: "#d7301f",
                8: "#b30000",
                9: "#7f0000"},

            "seqYellowOrangeRed": {
                1: "#ffffcc",
                2: "#ffeda0",
                3: "#fed976",
                4: "#feb24c",
                5: "#fd8d3c",
                6: "#fc4e2a",
                7: "#e31a1c",
                8: "#bd0026",
                9: "#800026"},

            "seqYellowOrangeBrown": {
                1: "#ffffe5",
                2: "#fff7bc",
                3: "#fee391",
                4: "#fec44f",
                5: "#fe9929",
                6: "#ec7014",
                7: "#cc4c02",
                8: "#993404",
                9: "#662506"},

             "seqPurples": {
                1: "#fcfbfd",
                2: "#efedf5",
                3: "#dadaeb",
                4: "#bcbddc",
                5: "#9e9ac8",
                6: "#807dba",
                7: "#6a51a3",
                8: "#54278f",
                9: "#3f007d"},

             "seqBlues": {
                1: "#f7fbff",
                2: "#deebf7",
                3: "#c6dbef",
                4: "#9ecae1",
                5: "#6baed6",
                6: "#4292c6",
                7: "#2171b5",
                8: "#08519c",
                9: "#08306b"},

             "seqGreens": {
                1: "#f7fcf5",
                2: "#e5f5e0",
                3: "#c7e9c0",
                4: "#a1d99b",
                5: "#74c476",
                6: "#41ab5d",
                7: "#238b45",
                8: "#006d2c",
                9: "#00441b"},

            "seqOranges": {
                1: "#fff5eb",
                2: "#fee6ce",
                3: "#fdd0a2",
                4: "#fdae6b",
                5: "#fd8d3c",
                6: "#f16913",
                7: "#d94801",
                8: "#a63603",
                9: "#7f2704"},

            "seqReds": {
                1: "#fff5f0",
                2: "#fee0d2",
                3: "#fcbba1",
                4: "#fc9272",
                5: "#fb6a4a",
                6: "#ef3b2c",
                7: "#cb181d",
                8: "#a50f15",
                9: "#67000d"},

            "seqGreys": {
                1: "#ffffff",
                2: "#f0f0f0",
                3: "#d9d9d9",
                4: "#bdbdbd",
                5: "#969696",
                6: "#737373",
                7: "#525252",
                8: "#252525",
                9: "#000000"},

            "divOrangePurple": {
                1: "#7f3b08",
                2: "#b35806",
                3: "#e08214",
                4: "#fdb863",
                5: "#fee0b6",
                6: "#f7f7f7",
                7: "#d8daeb",
                8: "#b2abd2",
                9: "#8073ac",
                10: "#542788",
                11: "#2d004b"},

            "divBrownBlueGreen": {
                1: "#543005",
                2: "#8c510a",
                3: "#bf812d",
                4: "#dfc27d",
                5: "#f6e8c3",
                6: "#f5f5f5",
                7: "#c7eae5",
                8: "#80cdc1",
                9: "#35978f",
                10: "#01665e",
                11: "#003c30"},

            "divPurpleGreen": {
                1: "#40004b",
                2: "#762a83",
                3: "#9970ab",
                4: "#c2a5cf",
                5: "#e7d4e8",
                6: "#f7f7f7",
                7: "#d9f0d3",
                8: "#a6dba0",
                9: "#5aae61",
                10: "#1b7837",
                11: "#00441b"},

            "divPinkGreen": {
                1: "#8e0152",
                2: "#c51b7d",
                3: "#de77ae",
                4: "#f1b6da",
                5: "#fde0ef",
                6: "#f7f7f7",
                7: "#e6f5d0",
                8: "#b8e186",
                9: "#7fbc41",
                10: "#4d9221",
                11: "#276419"},

            "divRedBlue": {
                1: "#67001f",
                2: "#b2182b",
                3: "#d6604d",
                4: "#f4a582",
                5: "#fddbc7",
                6: "#f7f7f7",
                7: "#d1e5f0",
                8: "#92c5de",
                9: "#4393c3",
                10: "#2166ac",
                11: "#053061"},

            "divRedGrey": {
                1: "#67001f",
                2: "#b2182b",
                3: "#d6604d",
                4: "#f4a582",
                5: "#fddbc7",
                6: "#ffffff",
                7: "#e0e0e0",
                8: "#bababa",
                9: "#878787",
                10: "#4d4d4d",
                11: "#1a1a1a"},

            "divRedYellowBlue": {
                1: "#a50026",
                2: "#d73027",
                3: "#f46d43",
                4: "#fdae61",
                5: "#fee090",
                6: "#ffffbf",
                7: "#e0f3f8",
                8: "#abd9e9",
                9: "#74add1",
                10: "#4575b4",
                11: "#313695"},

            "divRedYellowGreen": {
                1: "#a50026",
                2: "#d73027",
                3: "#f46d43",
                4: "#fdae61",
                5: "#fee08b",
                6: "#ffffbf",
                7: "#d9ef8b",
                8: "#a6d96a",
                9: "#66bd63",
                10: "#1a9850",
                11: "#006837"},

            "divSpectral": {
                1: "#9e0142",
                2: "#d53e4f",
                3: "#f46d43",
                4: "#fdae61",
                5: "#fee08b",
                6: "#ffffbf",
                7: "#e6f598",
                8: "#abdda4",
                9: "#66c2a5",
                10: "#3288bd",
                11: "#5e4fa2"},

            "qualSet1": {
                1: "#e41a1c",
                2: "#377eb8",
                3: "#4daf4a",
                4: "#984ea3",
                5: "#ff7f00",
                6: "#ffff33",
                7: "#a65628",
                8: "#f781bf",
                9: "#999999"},

            "qualSet2": {
                1: "#66c2a5",
                2: "#fc8d62",
                3: "#8da0cb",
                4: "#e78ac3",
                5: "#a6d854",
                6: "#ffd92f",
                7: "#e5c494",
                8: "#b3b3b3"},

            "qualSet3": {
                1: "#8dd3c7",
                2: "#ffffb3",
                3: "#bebada",
                4: "#fb8072",
                5: "#80b1d3",
                6: "#fdb462",
                7: "#b3de69",
                8: "#fccde5",
                9: "#d9d9d9",
                10: "#bc80bd",
                11: "#ccebc5",
                12: "#ffed6f"},

            "qualPastel1": {
                1: "#fbb4ae",
                2: "#b3cde3",
                3: "#ccebc5",
                4: "#decbe4",
                5: "#fed9a6",
                6: "#ffffcc",
                7: "#e5d8bd",
                8: "#fddaec",
                9: "#f2f2f2"},

            "qualPastel2": {
                1: "#b3e2cd",
                2: "#fdcdac",
                3: "#cbd5e8",
                4: "#f4cae4",
                5: "#e6f5c9",
                6: "#fff2ae",
                7: "#f1e2cc",
                8: "#cccccc"},

            "qualDark": {
                1: "#1b9e77",
                2: "#d95f02",
                3: "#7570b3",
                4: "#e7298a",
                5: "#66a61e",
                6: "#e6ab02",
                7: "#a6761d",
                8: "#666666"},

            "qualPaired": {
                1: "#a6cee3",
                2: "#1f78b4",
                3: "#b2df8a",
                4: "#33a02c",
                5: "#fb9a99",
                6: "#e31a1c",
                7: "#fdbf6f",
                8: "#ff7f00",
                9: "#cab2d6",
                10: "#6a3d9a",
                11: "#ffff99",
                12: "#b15928"},

            "qualAccent": {
                1: "#7fc97f",
                2: "#beaed4",
                3: "#fdc086",
                4: "#ffff99",
                5: "#386cb0",
                6: "#f0027f",
                7: "#bf5b17",
                8: "#666666"}
        }

    def _demo(self):
        """Draw all color maps."""
        import matplotlib.pyplot as plt
        map_names = sorted(self.maps.keys())
        num_maps = len(map_names)
        plot_num = 1
        num_cols = 13
        Xs = []
        Ys = []
        cols = []
        y = -1
        for map_name in map_names:
            y += 1
            for col in range(1, num_cols):
                try:
                    colour = self.maps[map_name][col]
                    Xs.append(col - 1)
                    Ys.append(y)
                    cols.append(colour)
                except KeyError:
                    pass
                plot_num += 1

        fig = plt.figure(figsize=(num_cols / 2, num_maps / 2.5), facecolor='w')
        ax = fig.add_subplot(111)
        ax.scatter(Xs,
                   Ys,
                   c=cols,
                   s=400,
                   edgecolors='none',
                   marker='h'
                   )
        plt.ylim((-1, num_maps)),
        plt.xticks(list(range(min(Xs), max(Xs) + 1)))
        plt.yticks(list(range(min(Ys), max(Ys) + 1)))
        ax.set_yticklabels(map_names)
        ax.set_xticklabels(list(range(1, num_cols)))
        plt.tight_layout()
        plt.savefig("cb2.png", dpi=300, format='png')
        plt.close(fig)
        del fig

if __name__ == "__main__":
    cb = ColorBrewer()
    cb._demo()
