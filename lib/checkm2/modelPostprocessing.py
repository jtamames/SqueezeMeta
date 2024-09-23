from checkm2.defaultValues import DefaultValues
import scipy
from scipy import sparse
from scipy.sparse import csr_matrix
import numpy as np
import logging
import sys
import os
import pickle
import pandas as pd


class modelProcessor:

    def __init__(self, threads=1):
        # self.logger = logging.getLogger('timestamp')

        try:
            self.ref_data = scipy.sparse.load_npz(DefaultValues.REF_DATA_LOCATION)
        except Exception as e:
            logging.error("Error: Reference data could not be loaded: {}".format(e))
            sys.exit(1)

        self.threads = threads
        self.reduced_cutoff = DefaultValues.AA_RATIO_COMPLETENESS_CUTOFF

        # Returns cosine similarity matrix

    def __calculate_sparse_CSM(self, A, B):

        # Reshape is require to deal with scipy's deviation from numpy calculations
        num = np.dot(A, B.T)

        p1 = np.sqrt(np.sum(A.power(2), axis=1))[:, np.newaxis]
        p1 = p1.reshape(p1.shape[0], p1.shape[1])

        p2 = np.sqrt(np.sum(B.power(2), axis=1))[np.newaxis, :]
        p2 = p2.reshape(p2.shape[0], p2.shape[1])

        return np.array( (num / (p1 * p2)).toarray() ) # FPS work with greater scipy version

    def __calculate_cosine_similarity(self, feature_vector):
        # Todo - if input is big, we might need to chunk it so as not to overload RAM (depends on ref data size). 

        csm = self.__calculate_sparse_CSM(self.ref_data[:, :20021], csr_matrix(feature_vector[:, :20021]))
        # return array of closest matches in ref database

        return np.amax(csm, axis=0)

    def cosine_decider(self, row):
        
        general = row['General']
        specific = row['Specific']
        cosine = row['Cosine_Similarity']
        AA_ratio = row['AA_Ratio']
        novelty_ratio = general / (cosine ** 2)

        meancomp = (general + specific) / 2


        if (meancomp < 55) and (AA_ratio < self.reduced_cutoff):
            #unfamiliar highly reduced genome - therefore needs to use general model
            return general, 'Gradient Boost (General Model)'
            
        else:
          if meancomp > 90:
            if novelty_ratio < 160:
                return specific, 'Neural Network (Specific Model)'
            else:
                return general, 'Gradient Boost (General Model)'
           
          if meancomp > 80:
              if novelty_ratio < 165:
                  return specific, 'Neural Network (Specific Model)'
              else:
                  return general, 'Gradient Boost (General Model)'
          if meancomp > 70:
              if novelty_ratio < 165:
                  return specific, 'Neural Network (Specific Model)'
              else:
                  return general, 'Gradient Boost (General Model)'
          if meancomp > 60:
              if novelty_ratio < 170:
                  return specific, 'Neural Network (Specific Model)'
              else:
                  return general, 'Gradient Boost (General Model)'
          if meancomp > 50:
              if novelty_ratio < 175:
                  return specific, 'Neural Network (Specific Model)'
              else:
                  return general, 'Gradient Boost (General Model)'
          if meancomp > 40:
              if novelty_ratio < 175:
                  return specific, 'Neural Network (Specific Model)'
              else:
                  return general, 'Gradient Boost (General Model)'
          else:
              return specific, 'Neural Network (Specific Model)'
        

    def calculate_general_specific_ratio(self, AA_counts, feature_vector, general_comp, contamination, specific_comp):



        csm_array = self.__calculate_cosine_similarity(csr_matrix(feature_vector))

        mean_completeness = (general_comp + specific_comp) / 2
        completeness_AA_ratio = AA_counts / mean_completeness

        comp_results = pd.DataFrame({'General': general_comp, 
                                     'Specific': specific_comp, 
                                     'Cosine_Similarity': csm_array, 
                                     'AA_Ratio': completeness_AA_ratio})
                                     

        comp_results[['CheckM2_Completeness', 'Model_Chosen']] = comp_results.apply(self.cosine_decider, axis=1, result_type='expand')
        
        return comp_results['CheckM2_Completeness'].values, contamination, comp_results['Model_Chosen'].values, comp_results['Cosine_Similarity'].values

        # specific_bool_mask = (csm_array > DefaultValues.COSINE_SIMILARITY_SPECIFIC_MODEL_CUTOFF)
    
        # final_completeness = np.where(specific_bool_mask, specific_comp, general_comp)
        # final_contamination = np.where(specific_bool_mask, specific_cont, general_cont)
    
        # models_chosen = ['Specific (NeuralNetwork)' if specific else 'General (GradientBoost)' for specific in specific_bool_mask]
    
        #return final_completeness, final_contamination, models_chosen, csm_array

