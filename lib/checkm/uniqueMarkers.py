###############################################################################
#
# uniqueMarkers.py - find a set of markers that are descriptive for a taxonomy
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

import argparse
import sqlite3
import re
from collections import defaultdict


def parseTaxonomy(taxonomy):
    tax = re.compile('[;,:]?\s?[kpcofgs]__|[;,:]')
    return tax.split(taxonomy)[1:]


def getTaxonId(cursor, *args):
    ranks = ['Domain', 'Phylum', 'Class', '"Order"', 'Family', 'Genus', 'Species']
    query = []
    for rank, value in zip(ranks, args):
        query.append(' %s = \'%s\' ' % (rank, value))
    query.append(' %s IS NULL' % ranks[len(args)])
    query_string = 'AND'.join(query)
    result = cursor.execute('SELECT Id, "Count" FROM taxons WHERE %s' % query_string)
    return result.fetchall()


def getOppositeRankSpecificTaxonId(cursor, *args):
    ''' Get all other taxon lineages at the same level as the requested taxon
    '''
    ranks = ['Domain', 'Phylum', 'Class', '"Order"', 'Family', 'Genus', 'Species']
    query = []
    for rank, value in zip(ranks, args[:-1]):
        query.append(' %s = \'%s\' ' % (rank, value))
    query.append(' %s != \'%s\' ' % (ranks[len(args) - 1], args[-1]))
    query.append(' %s IS NULL' % ranks[len(args)])
    query_string = 'AND'.join(query)
    print query_string
    result = cursor.execute('SELECT Id, "Count" FROM taxons WHERE %s' % query_string)
    return result.fetchall()


def getOppositeRankInspecificTaxonId(cursor, *args):
    ''' Get all other taxon lineages at the same level as the requested taxon
    '''
    ranks = ['Domain', 'Phylum', 'Class', '"Order"', 'Family', 'Genus', 'Species']
    query = []
    for rank, value in zip(ranks, args):
        query.append(' %s != \'%s\' ' % (rank, value))
    # query.append(' %s IS NULL' % ranks[len(args)])
    query_string = query[-1]
    result = cursor.execute('SELECT Id, "Count" FROM taxons WHERE %s' % query_string)
    return result.fetchall()


def getMarkersFromTaxon(cursor, taxid):
    result = cursor.execute('''SELECT Marker, "Count" FROM marker_mapping WHERE Taxon = ?''', (taxid,))
    return result.fetchall()


def getMarkersNotInTaxon(cursor, taxid):
    result = cursor.execute('''SELECT Marker, "Count" FROM marker_mapping WHERE Taxon != ?''', (taxid,))
    return result.fetchall()


def countAllGenomes(cursor):
    result = cursor.execute('''SELECT Id, "Count" FROM taxons''')
    return result.fetchall()


def countGenomesInTaxon(cursor, taxId):
    result = cursor.execute('''SELECT "Count" FROM taxons WHERE Id = ?''', (taxId,))
    return result.fetchone()[0]


def getDescriptiveMarkers(cursor, markers):
    result = cursor.execute('''SELECT Acc, Name FROM markers WHERE Id = ?''', markers)
    return result.fetchone()


def doWork(args):
    taxon_ranks = parseTaxonomy(args.taxonomy)
    con = sqlite3.connect(args.database)
    cur = con.cursor()

    taxon_ids = getTaxonId(cur, *taxon_ranks)
    if len(taxon_ids) > 1:
        raise RuntimeError("Taxon string returns more than one lineage "\
                "please be more specific")
    else:
        tax_id, tax_count = taxon_ids[0]
        all_markers = getMarkersFromTaxon(cur, tax_id)
        marker_in_taxon_mapping = {}
        for (Id, count) in all_markers:
            if float(count) / float(tax_count) >= args.include:
                marker_in_taxon_mapping[Id] = float(count) / float(tax_count)

        opposite_taxons = getOppositeRankInspecificTaxonId(cur, *taxon_ranks)
        markers_from_others = defaultdict(int)
        others_total_count = 0
        for (Id, count) in opposite_taxons:
            others_total_count += count
            this_taxon_count = getMarkersFromTaxon(cur, Id)
            for (Id, count) in this_taxon_count:
                markers_from_others[Id] += count

        descriptive_markers = []
        for marker_id, _ in marker_in_taxon_mapping.items():
            if marker_id in markers_from_others:
                fraction_in_others = float(markers_from_others[marker_id]) / float(others_total_count)
                if fraction_in_others <= args.exclude:
                    descriptive_markers.append((marker_id,))
            else:
                # not found in anything else
                descriptive_markers.append((marker_id,))

        des_markers = []
        for i in descriptive_markers:
            des_markers.append(getDescriptiveMarkers(cur, i))

        for des_acc, des_name in des_markers:
            print des_acc, des_name

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-t', '--taxon-string',
            default='k__Bacteria;p__Proteobacteria',
            dest='taxonomy', help="specify the taxonomy")
    parser.add_argument('-i', '--inclusive-percent', dest='include',
            default=0.95, help="The required percentage of member of the "\
            "specified taxonomy must have a particular marker")
    parser.add_argument('-e', '--exclusive-percent', dest='exclude',
            default=0.05, help="The maximum prevelence of the marker in "\
            "all other non-target taxons")
    parser.add_argument('-d', '--database', dest='database',
        default='markers.db', help='specify path to database')

    # parse the arguments
    args = parser.parse_args()

    # do what we came here to do
    doWork(args)
