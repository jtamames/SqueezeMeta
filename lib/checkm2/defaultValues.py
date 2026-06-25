import os
from pathlib import Path


class DefaultValues():
    """Default values for filenames and common constants."""
    DATA_PATH =  os.path.join(os.path.dirname(__file__), 'data')
    MODEL_PATH = os.path.join(os.path.dirname(__file__), 'models')
    VERSION_PATH = os.path.join(os.path.dirname(__file__), 'version')
    TESTRUN_GENOMES = os.path.join(os.path.dirname(__file__), 'testrun')

    PRODIGAL_FOLDER_NAME = 'protein_files'

    #LOGNAME = 'log.txt'

    DIAMOND_QUERY_COVER = 80
    DIAMOND_SUBJECT_COVER = 80
    DIAMOND_PERCENT_ID = 30
    DIAMOND_EVALUE = 1e-05

    DIAMOND_DEFAULT_CHUNK_SIZE = 250

    KO_FEATURE_VECTOR_CHUNK = 250

    DIAMOND_HEADER_SEPARATOR = 'Ω'
    
    AA_RATIO_COMPLETENESS_CUTOFF = 1500

    MODEL_DIVERGENCE_WARNING_THRESHOLD = 25

    DB_LOCATION_DEFINITION = os.path.join(VERSION_PATH, 'diamond_path.json')
    DB_VAR = "CHECKM2DB"
    try:
        DEFAULT_DB_INSTALL_LOCATION = os.path.join(str(Path.home()), 'databases')
    except:
        pass


    FEATURE_ORDER_LOCATION = os.path.join(DATA_PATH, 'feature_ordering.json')
    PATH_CATEGORY_MAPPING_LOCATION = os.path.join(DATA_PATH, 'kegg_path_category_mapping.json')
    MODULE_DEFINITION_LOCATION = os.path.join(DATA_PATH, 'module_definitions.json')

    GENERAL_MODEL_COMP_LOCATION = os.path.join(MODEL_PATH, 'general_model_COMP.gbm')
    SPECIFIC_MODEL_COMP_LOCATION = os.path.join(MODEL_PATH, 'specific_model_COMP.keras')

    MODEL_CONT_LOCATION = os.path.join(MODEL_PATH, 'model_CONT.gbm')
#    SPECIFIC_MODEL_CONT_LOCATION = os.path.join(MODEL_PATH, 'specific_model_CONT.hd5')

    SCALER_FILE_LOCATION = os.path.join(MODEL_PATH, 'scaler.sav')
    COSINE_TABLE_LOCATION = os.path.join(MODEL_PATH, 'cosine_table.pkl')
    REF_DATA_LOCATION = os.path.join(DATA_PATH, 'min_ref_rsdata_v1.npz')

    EXTERNAL_FILES_TO_VERIFY = [FEATURE_ORDER_LOCATION, PATH_CATEGORY_MAPPING_LOCATION,
                                MODULE_DEFINITION_LOCATION, GENERAL_MODEL_COMP_LOCATION,
                                SPECIFIC_MODEL_COMP_LOCATION, MODEL_CONT_LOCATION,
                                SCALER_FILE_LOCATION, REF_DATA_LOCATION, COSINE_TABLE_LOCATION]

