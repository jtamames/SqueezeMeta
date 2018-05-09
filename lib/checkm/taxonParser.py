###############################################################################
#
# taxonParser.py - parse taxonomic-specific marker sets
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

import logging
from collections import defaultdict

import checkm.prettytable as prettytable

from checkm.markerSets import BinMarkerSets, MarkerSet
from checkm.util.taxonomyUtils import taxonomicRanks, ranksByLevel, ranksByLabel

from checkm.defaultValues import DefaultValues


class TaxonParser():
    """Parse taxonomic-specific marker sets."""
    def __init__(self):
        self.logger = logging.getLogger()

    def readMarkerSets(self):
        taxonMarkerSets = defaultdict(dict)
        for line in open(DefaultValues.TAXON_MARKER_SETS):
            lineSplit = line.split('\t')
            rank = lineSplit[0]
            taxon = lineSplit[1]
            lineage = lineSplit[2]
            numGenomes = int(lineSplit[3])
            markerSet = eval(lineSplit[6].rstrip())
            
            ms = MarkerSet(ranksByLabel[rank], lineage, numGenomes, markerSet)
            ms.removeMarkers(DefaultValues.MARKERS_TO_EXCLUDE)

            taxonMarkerSets[rank][taxon] = ms

        return taxonMarkerSets

    def list(self, rankFilter='ALL'):
        """ List all available marker sets from the specified rank."""

        taxonMarkerSets = self.readMarkerSets()

        header = ['Rank', 'Taxon', '# genomes', '# marker genes', '# marker sets']
        pTable = prettytable.PrettyTable(header)
        pTable.align = 'c'
        pTable.align['Rank'] = 'l'
        pTable.align['Taxon'] = 'l'
        pTable.hrules = prettytable.FRAME
        pTable.vrules = prettytable.NONE

        for rank in taxonomicRanks:
            if rankFilter == 'ALL' or rankFilter == rank:
                for taxon in sorted(taxonMarkerSets[rank]):
                    markerSet = taxonMarkerSets[rank][taxon]

                    numMarkers, numMarkerSets = markerSet.size()
                    pTable.add_row([rank, taxon, markerSet.numGenomes, numMarkers, numMarkerSets])

        print ''
        print pTable.get_string()

    def markerSet(self, rank, taxon, markerFile):
        """Obtain specified taxonomic-specific marker set."""

        taxonMarkerSets = self.readMarkerSets()

        if rank not in taxonMarkerSets:
            self.logger.error('  Unrecognized taxonomic rank: ' + rank)
            return False
        elif taxon not in taxonMarkerSets[rank]:
            self.logger.error('  Unrecognized taxon: %s (in rank %s): ' % (taxon, rank))
            return False

        markerSet = taxonMarkerSets[rank][taxon]

        taxonomy = markerSet.lineageStr.split(';')[::-1]
        binMarkerSets = BinMarkerSets(taxon, BinMarkerSets.TAXONOMIC_MARKER_SET)
        for i, taxon in enumerate(taxonomy):
            if rank != 'life':
                rank = ranksByLevel[len(taxonomy) - i - 1]

            if rank == 'species':
                taxon = taxonomy[1] + ' ' + taxonomy[0]

            markerSet = taxonMarkerSets[rank][taxon]
            numMarkers, numMarkerSets = markerSet.size()
            self.logger.info('  Marker set for %s contains %d marker genes arranged in %d sets.' % (taxon, numMarkers, numMarkerSets))
            self.logger.info('    Marker set inferred from %d reference genomes.' % markerSet.numGenomes)

            markerSet.lineageStr = taxon
            binMarkerSets.addMarkerSet(markerSet)

        fout = open(markerFile, 'w')
        fout.write(DefaultValues.TAXON_MARKER_FILE_HEADER + '\n')
        binMarkerSets.write(fout)
        fout.close()

        return True
