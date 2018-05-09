#!/usr/bin/env python

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

"""
Useful methods for processing taxonomy strings.
"""

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2013'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '1.0.0'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

taxonomicRanks = ['life', 'domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']
taxonomicPrefixes = ['l__', 'k__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']

ranksByLabel = {'life':-1, 'domain': 0, 'phylum': 1, 'class': 2, 'order': 3, 'family': 4, 'genus': 5, 'species': 6, 'gg_id': 7}
ranksByLevel = {-1: 'life', 0: 'domain', 1: 'phylum', 2: 'class', 3: 'order', 4: 'family', 5: 'genus', 6: 'species', 7: 'gg_id'}
rankPrefixes = {-1: 'l__', 0: 'k__', 1: 'p__', 2: 'c__', 3: 'o__', 4: 'f__', 5: 'g__', 6: 's__', 7: 'id__'}


def appendTaxonomyRanks(taxonomy, ranks=7):
    """Append taxonomy prefix to list of taxa."""
    t = []
    for i in xrange(0, ranks):
        t.append(rankPrefixes[i] + taxonomy[i])

    return t


def readTaxonomy(taxonomyFile):
    """Read tab-seperated taxonomy file."""
    taxonIdToTaxonomy = {}
    for line in open(taxonomyFile):
        lineSplit = line.split('\t')
        taxonIdToTaxonomy[lineSplit[0]] = lineSplit[1].rstrip()

    return taxonIdToTaxonomy


def __parseTaxon(taxon):
    if '(' in taxon:
        taxonSplit = taxon.split('(')
        taxonId = taxonSplit[0]
        taxonId = taxonId.strip()
        bootstrapSupport = int(taxonSplit[1][0:taxonSplit[1].find(')')])
    else:
        taxonId = taxon.strip()
        bootstrapSupport = 0

    return taxonId, bootstrapSupport


def LCA(taxonomy1, taxonomy2):
    """Find lowest-common ancestor between two taxa lists."""
    lca = []
    for i in xrange(0, min(len(taxonomy1), len(taxonomy2))):
        t1, b1 = __parseTaxon(taxonomy1[i])
        t2, b2 = __parseTaxon(taxonomy2[i])

        if t1 != t2:
            if 'unmapped' in t1 or 'unmapped' in t2:
                lca.append(rankPrefixes[i] + 'unmapped')
            else:
                lca.append(rankPrefixes[i] + 'unclassified')
        else:
            if b1 == 0 and b2 == 0:
                lca.append(t1)
            else:
                lca.append(t1 + '(' + str(min(b1, b2)) + ')')

    return lca
