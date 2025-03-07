from checkm2 import fileManager
from checkm2.defaultValues import DefaultValues

import json
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'


class KeggCalculator:

    ''' Dependent on JSON-defined gene, pathway, category and module information

        Feature_order: Dictionary with ordered keys (categories) containing ordered values (feature columns)
        Path_category_mapping: pd Dataframe with Kegg_ID, Pathway, Category as columns
        module_definitions: Dictionary with module names as keys and corresponding KO's as values

        Features are ordered as follows:
       (['Metadata', 'KO_Genes', 'KO_Pathways', 'KO_Modules', 'KO_Categories'])
    '''
    def __init__(self):

        with open(DefaultValues.FEATURE_ORDER_LOCATION, 'r') as fo:
            self.feature_order = json.load(fo)

        with open(DefaultValues.PATH_CATEGORY_MAPPING_LOCATION, 'r') as pcm:
            self.path_category_mapping = pd.DataFrame(json.load(pcm))

        with open(DefaultValues.MODULE_DEFINITION_LOCATION, 'r') as md:
            self.module_definitions = json.load(md)


    def return_default_values_from_category(self, ID):
        return dict.fromkeys(self.feature_order[ID], 0)

    def return_proper_order(self, category):
        return self.feature_order[category]

    def calculate_KO_group(self, group, KO_gene_data):

        #only interested in presence/absence; also exclude names (last column)
        presence_absence_subset = KO_gene_data.iloc[:, :-1].copy()
        presence_absence_subset[presence_absence_subset > 1] = 1

        #get ordered sequence for each feature vector in group
        ordered_entries = self.return_default_values_from_category(group)
        feature_vectors = list(ordered_entries.keys())

        for vector in feature_vectors:

            #Go through category, determine total Kegg_IDs belonging to it
            length = len(self.path_category_mapping[self.path_category_mapping[group] == vector]['Kegg_ID'])

            #Count how many are present
            KO_gene_data[vector] = (presence_absence_subset[self.path_category_mapping[self.path_category_mapping[group] == vector]['Kegg_ID'].values].sum(axis=1)) / length

        #return only that part of dataframe we're interested in
        return KO_gene_data.iloc[:, -len(feature_vectors):]

    '''Module calculations differ from others as one gene can be part of multiple modules'''


    def calculate_module_completeness(self, KO_gene_data):
        presence_absence_subset = KO_gene_data.iloc[:, :-1].copy()
        presence_absence_subset[presence_absence_subset > 1] = 1

        modules = list(self.module_definitions.keys())

        for module in modules:
            length = len(self.module_definitions[module])

            #make sure all components of the module definition are accounted for
            self.module_definitions[module] = [x for x in self.module_definitions[module] if x in KO_gene_data.columns]

            KO_gene_data[module] = (KO_gene_data[self.module_definitions[module]].sum(axis=1)) / length

        return KO_gene_data.iloc[:, -len(modules):]



