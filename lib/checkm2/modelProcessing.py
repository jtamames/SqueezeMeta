import warnings

warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=DeprecationWarning)

import pandas as pd
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)

from checkm2.defaultValues import DefaultValues
# import xgboost as xgb
import lightgbm as lgb
import os

#make sure we're only using CPUs as GPUs can throw weird errors and is not worth the minor speed advantage
os.environ['CUDA_VISIBLE_DEVICES'] = '-1'

#from tensorflow import keras
import keras
from keras.models import load_model
os.environ["KERAS_BACKEND"] = "tensorflow"

from sklearn.preprocessing import MinMaxScaler
import pickle
import logging
import sys
import os


class modelProcessor:

    def __init__(self, threads):



        self.nthreads = threads

        try:
            self.general_model_comp = lgb.Booster(model_file=DefaultValues.GENERAL_MODEL_COMP_LOCATION)
            self.model_cont = lgb.Booster(model_file=DefaultValues.MODEL_CONT_LOCATION)

            self.specific_model_comp_nn = load_model(DefaultValues.SPECIFIC_MODEL_COMP_LOCATION)

            self.minmax_scaler = pickle.load(open(DefaultValues.SCALER_FILE_LOCATION, 'rb'))
            

            if logging.root.level == logging.DEBUG:
                self.verbosity = 1
            else:
                self.verbosity = 0

        except Exception as e:
            logging.error("Saved models could not be loaded: {}".format(e))
            sys.exit(1)

    def run_prediction_general(self, vector_array):
        
        #TODO: make sure runs on 1 sample

        comp_predictions = self.general_model_comp.predict(vector_array, n_jobs=self.nthreads)
        comp_predictions[comp_predictions > 100] = 100

        cont_predictions = self.model_cont.predict(vector_array, n_jobs=self.nthreads)

        comp_predictions[comp_predictions < 0] = 0
        cont_predictions[cont_predictions < 0] = 0

        return comp_predictions.flatten(), cont_predictions.flatten()

    def run_prediction_specific(self, vector_array, specific_model_vector_len):

        scaled_vector = self.minmax_scaler.transform(vector_array)

        # re-shape into keras-cnn-appropriate array
        scaled_vector = scaled_vector.reshape(scaled_vector.shape[0], scaled_vector.shape[1], 1)

        # only using genes for specific predictions

        comp_predictions = self.specific_model_comp_nn.predict(scaled_vector[:, :specific_model_vector_len],
                                                               verbose=self.verbosity)

        # as we're using sigmoid output for completeness model, convert to 100-scale
        comp_predictions = comp_predictions * 100


        comp_predictions[comp_predictions < 0] = 0

        return comp_predictions.flatten(), scaled_vector.reshape(scaled_vector.shape[0], scaled_vector.shape[1])

