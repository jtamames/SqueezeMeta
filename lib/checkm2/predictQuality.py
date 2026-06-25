from checkm2 import modelProcessing
from checkm2 import metadata
from checkm2 import prodigal
from checkm2 import diamond
from checkm2.defaultValues import DefaultValues
from checkm2.versionControl import VersionControl
from checkm2 import keggData
from checkm2 import modelPostprocessing
from checkm2 import fileManager

import os
import multiprocessing as mp
import numpy as np
import shutil
import sys
import logging
import pandas as pd
import tarfile

# For unnessesary tensorflow warnings:
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
logging.getLogger('tensorflow').setLevel(logging.FATAL)



class Predictor():
    def __init__(self, bin_folder, outdir, bin_extension='.fna', threads=1, lowmem=False, tempDBloc=None):

        self.bin_folder = bin_folder
        self.bin_extension = bin_extension
        self.bin_files = self.__setup_bins()

        self.output_folder = outdir
        self.prodigal_folder = os.path.join(self.output_folder, DefaultValues.PRODIGAL_FOLDER_NAME)
        fileManager.make_sure_path_exists(self.prodigal_folder)
        
        self.lowmem = lowmem
        if self.lowmem:
          logging.info('Running in low-memory mode.')


        self.total_threads = threads

        logging.debug('Verifying internal checksums for all models, scalers and reference data.')
        #if VersionControl().checksum_version_validate() is False:
        #    logging.error('Could not verify internal model checksums. Please re-download CheckM2.')
        #    sys.exit(1)

        logging.debug('Verifying DIAMOND DB installation path.')

        if tempDBloc is not None:
            self.diamond_path = tempDBloc
        else:
            self.diamond_path = fileManager.DiamondDB().get_DB_location()

        if self.diamond_path == None or self.diamond_path == '' or self.diamond_path == 'Not Set':
            logging.error("Please download and install the CheckM2 database first (see 'checkm2 database -h')")
            sys.exit(1)

        fileManager.check_if_file_exists(self.diamond_path)


    def __setup_bins(self):
        bin_files = []
        if self.bin_folder is not None:
            all_files = os.listdir(self.bin_folder)
            for f in all_files:
                if f.endswith(self.bin_extension):
                    binFile = os.path.join(self.bin_folder, f)
                    if os.stat(binFile).st_size == 0:
                        logging.warning("Skipping bin {} as it has a size of 0 bytes.".format(f))
                    elif tarfile.is_tarfile(binFile):
                      logging.warning('Skipping bin {} as tar archives are not supported.'.format(binFile))
                    else:
                        bin_files.append(binFile)

        if not bin_files:
            logging.error("No bins found. Check the extension (-x) used to identify bins.")
            sys.exit(1)

        return sorted(bin_files)

    def prediction_wf(self, genes_supplied=False, mode='auto', debug_cos=False,
                      dumpvectors=False, stdout=False, resume=False, remove_intermediates=False, ttable=None):

        #make sure models can be loaded without problems
        modelProc = modelProcessing.modelProcessor(self.total_threads)

        #make sure diamond is set up and ready to go
        diamond_search = diamond.DiamondRunner(self.total_threads, self.output_folder, self.lowmem, self.diamond_path)



        ''' 1: Call genes and automatically determine coding table'''
        if resume:
            logging.info('Re-using protein files from output directory: {}'.format(self.prodigal_folder,))
            prodigal_files = [os.path.join(self.prodigal_folder, bin_file) for bin_file in os.listdir(self.prodigal_folder)]

        elif not genes_supplied:
            used_ttables, coding_density, \
            N50, avg_gene_len, \
            total_bases, cds_count, \
            GC, totalContigs, maxContigLen = self.__run_prodigal(ttable)

            prodigal_files, used_ttables = fileManager.verify_prodigal_output(self.prodigal_folder, used_ttables, self.bin_extension)

        else:
            logging.info('Using user-supplied protein files.')
            prodigal_files = []
            for bin in self.bin_files:
                shutil.copyfile(bin, os.path.join(self.prodigal_folder, os.path.splitext(os.path.basename(bin))[0]))
                prodigal_files.append(bin)

        ''' 2: Calculate genome metadata from protein files'''

        metadata_df = self.__calculate_metadata(prodigal_files)
        metadata_df = pd.concat(metadata_df.values())
        metadata_df.reset_index(drop=True, inplace=True)

        # make sure metadata is arranged correctly
        metadata_order = keggData.KeggCalculator().return_proper_order('Metadata')
        metadata_order.insert(0, 'Name')
        metadata_df = metadata_df[metadata_order]

        ''' 3: Determine all KEGG annotations of input genomes using DIAMOND blastp'''

        if resume:
            logging.info("Reusing DIAMOND output from output directory: {}".format(diamond_search.diamond_out))
                
            diamond_out = [x for x in os.listdir(diamond_search.diamond_out) if x.startswith('DIAMOND_RESULTS')]
            if len(diamond_out) == 0:
                logging.error("No DIAMOND outputs have been found in {}. Resuming is not possible.".format(diamond_search.diamond_out))
                exit(1)
        else:
            diamond_out = diamond_search.run(prodigal_files)

        ### MOVED

        logging.info('Processing DIAMOND output')
        # concatenate all results even if only one
        results = pd.concat([pd.read_csv(os.path.join(diamond_search.diamond_out, entry), sep='\t', usecols=[0, 1],
                                         names=['header', 'annotation']) for entry in diamond_out])

        if len(results) < 1:
            logging.error('No DIAMOND annotation was generated. Exiting')
            sys.exit(1)

        # Split columns into usable series
        results[['GenomeName', 'ProteinID']] = results['header'].str.split(diamond_search.separator, n=1, expand=True)
        results[['Ref100_hit', 'Kegg_annotation']] = results['annotation'].str.split('~', n=1, expand=True)



        ''' Get a list of default KO id's from data
            Available categories are the keys in DefaultValues.feature_ordering
            Here, returns an ordered set of KEGG ID's and sets to 0 
        '''
        KeggCalc = keggData.KeggCalculator()
        defaultKOs = KeggCalc.return_default_values_from_category('KO_Genes')

        # Remove from results any KOs we're not currently using
        results = results[results['Kegg_annotation'].isin(defaultKOs.keys())]

        # Update counts per genome

        full_name_list = metadata_df['Name'].values
        #kegg_genome_list = []

        annot_dict = dict(
            zip(sorted(results['GenomeName'].unique()), [x for _, x in results.groupby(results['GenomeName'])]))

        logging.info('Predicting completeness and contamination using ML models.')

        names, final_comps, final_conts, models_chosen, csm_arrays, \
        general_results_comp, specific_results_comp = [], [], [], [], [], [], []

        chunk_counter = 0

        for i in range(0, len(full_name_list), DefaultValues.KO_FEATURE_VECTOR_CHUNK):
            sublist = full_name_list[i:i + DefaultValues.KO_FEATURE_VECTOR_CHUNK]
            chunk_counter += 1


            parsed_diamond_results, ko_list_length = diamond_search.process_diamond_output(defaultKOs, annot_dict, sublist)


            parsed_diamond_results.sort_values(by='Name', inplace=True)

            sub_metadata = metadata_df[metadata_df['Name'].isin(sublist)]

            sub_metadata.sort_values(by='Name', inplace=True)
            parsed_diamond_results.sort_values(by='Name', inplace=True)
            parsed_diamond_results.reset_index(drop=True, inplace=True)
            sub_metadata.reset_index(drop=True, inplace=True)

            names.append(parsed_diamond_results['Name'].values)

            # delete duplicate 'name' column and merge
            del parsed_diamond_results['Name']

            feature_vectors = pd.concat([sub_metadata[sub_metadata['Name'].isin(sublist)], parsed_diamond_results], axis=1)
            #print(feature_vectors.shape)

            ''' 4: Call general model & specific models and derive predictions'''

            vector_array = feature_vectors.iloc[:, 1:].values.astype(float)


            general_result_comp, general_result_cont = modelProc.run_prediction_general(vector_array)

            specific_model_vector_len = (ko_list_length + len(
                metadata_order)) - 1  # -1 = without name TODO a bit ugly - maybe just calculate length on setup somewhere


            # also retrieve scaled data for CSM calculations
            specific_result_comp, scaled_features = modelProc.run_prediction_specific(vector_array, specific_model_vector_len)

            final_conts.append(general_result_cont)
            general_results_comp.append(general_result_comp)
            specific_results_comp.append(specific_result_comp)


            ''' 5: Determine any substantially complete genomes similar to reference genomes and fine-tune predictions'''

            if not mode == 'specific' or not mode == 'general':
                #logging.info('Using cosine simlarity to reference data to select appropriate predictor model.')

                postProcessor = modelPostprocessing.modelProcessor(self.total_threads)
                final_comp, final_cont, model_chosen, csm_array = postProcessor.calculate_general_specific_ratio(
                    vector_array[:, 20],
                    scaled_features,
                    general_result_comp,
                    general_result_cont,
                    specific_result_comp)

                final_comps.append(final_comp)
                models_chosen.append(model_chosen)
                csm_arrays.append(csm_array)

            if dumpvectors:
                dumpfile = os.path.join(self.output_folder, f'feature_vectors_{chunk_counter}.pkl')
                feature_vectors.to_pickle(dumpfile, protocol=4)

        logging.info('Parsing all results and constructing final output table.')


        #flatten lists
        names = [item for sublist in names for item in sublist]
        final_comps = [item for sublist in final_comps for item in sublist]
        final_conts = [item for sublist in final_conts for item in sublist]
        models_chosen = [item for sublist in models_chosen for item in sublist]
        csm_arrays = [item for sublist in csm_arrays for item in sublist]
        general_results_comp = [item for sublist in general_results_comp for item in sublist]
        specific_results_comp = [item for sublist in specific_results_comp for item in sublist]


        final_results = pd.DataFrame({'Name':names})

        if mode == 'both':
            final_results['Completeness_General'] = np.round(general_results_comp, 2)
            final_results['Contamination'] = np.round(final_conts, 2)
            final_results['Completeness_Specific'] = np.round(specific_results_comp, 2)
            final_results['Completeness_Model_Used'] = models_chosen

        elif mode == 'auto':
            final_results['Completeness'] = np.round(final_comps, 2)
            final_results['Contamination'] = np.round(final_conts, 2)
            final_results['Completeness_Model_Used'] = models_chosen

        elif mode == 'general':
            final_results['Completeness_General'] = np.round(general_results_comp, 2)
            final_results['Contamination'] = np.round(final_conts, 2)

        elif mode == 'specific':
            final_results['Completeness_Specific'] = np.round(specific_results_comp, 2)
            final_results['Contamination'] = np.round(final_conts, 2)

        else:
            logging.error('Programming error in model choice')
            sys.exit(1)

        if not genes_supplied and not resume:
            final_results['Translation_Table_Used'] = final_results['Name'].apply(lambda x: used_ttables[x])
            final_results['Coding_Density'] = final_results['Name'].apply(lambda x: np.round(coding_density[x], 3))
            final_results['Contig_N50'] = final_results['Name'].apply(lambda x: int(N50[x]))
            final_results['Average_Gene_Length'] = final_results['Name'].apply(lambda x: avg_gene_len[x])
            final_results['Genome_Size'] = final_results['Name'].apply(lambda x: total_bases[x])
            final_results['GC_Content'] = final_results['Name'].apply(lambda x: np.round(GC[x], 2))
            final_results['Total_Coding_Sequences'] = final_results['Name'].apply(lambda x: cds_count[x])
            final_results['Total_Contigs'] = final_results['Name'].apply(lambda x: totalContigs[x])
            final_results['Max_Contig_Length'] = final_results['Name'].apply(lambda x: maxContigLen[x])


        if debug_cos is True:
            final_results['Cosine_Similarity'] = np.round(csm_arrays, 2)


            
        #Flag any substantial divergences in completeness predictions
        additional_notes = self.__flag_divergent_predictions(general=general_results_comp, specific=specific_results_comp)
        
        final_results['Additional_Notes'] = additional_notes

        final_file = os.path.join(self.output_folder, 'quality_report.tsv')
        final_results.to_csv(final_file, sep='\t', index=False)
        
        if stdout:
            print(final_results.to_string(index=False, float_format=lambda x: '%.2f' % x))

        if remove_intermediates:
            shutil.rmtree(self.prodigal_folder)
            shutil.rmtree(diamond_search.diamond_out)

        logging.info('CheckM2 finished successfully.')

    def __flag_divergent_predictions(self, general, specific, threshold=DefaultValues.MODEL_DIVERGENCE_WARNING_THRESHOLD):
    
        compare = pd.DataFrame({'General':general, 'Specific':specific})
        compare['Difference'] = compare.apply(lambda row: abs(row['General'] - row['Specific']), axis=1)
        compare['Additional_Notes'] = compare.apply(lambda row: 'None' if row['Specific'] < 50 or row['Difference'] < threshold else \
        'Low confidence prediction - substantial ({}%) disagreement between completeness prediction models'.format(int(row['Difference'])), axis=1)
        
        return compare['Additional_Notes'].values

    def __set_up_prodigal_thread(self, queue_in, queue_out, ttable, used_ttable, coding_density,
                                 N50, avg_gene_len, total_bases, cds_count, GC, totalContigs, maxContigLen):

        while True:
            bin = queue_in.get(block=True, timeout=None)
            if bin == None:
                break

            prodigal_thread = prodigal.ProdigalRunner(self.prodigal_folder, bin)
            binname, selected_coding_table, c_density, \
            v_N50, v_avg_gene_len, v_total_bases, v_cds_count, \
            v_GC, v_totalContigs, v_maxContigLen = prodigal_thread.run(bin, ttable)

            used_ttable[binname] = selected_coding_table
            coding_density[binname] = c_density
            N50[binname] = v_N50
            avg_gene_len[binname] = v_avg_gene_len
            total_bases[binname] = v_total_bases
            GC[binname] = v_GC
            cds_count[binname] = v_cds_count
            totalContigs[binname] = v_totalContigs
            maxContigLen[binname] = v_maxContigLen

            queue_out.put((bin, selected_coding_table, coding_density, N50, avg_gene_len, total_bases, cds_count,
                           GC, totalContigs, maxContigLen))

    def __reportProgress(self, total_bins, queueIn):
        """Report number of processed bins."""

        processed = 0

        while True:
            bin, selected_coding_table, coding_density, N50, \
            avg_gene_len, total_bases, cds_count, GC, totalContigs, maxContigLen = queueIn.get(block=True, timeout=None)
            if bin == None:
                if logging.root.level == logging.INFO or logging.root.level == logging.DEBUG:
                    sys.stdout.write('\n')
                    sys.stdout.flush()
                break

            processed += 1

            if logging.root.level == logging.INFO or logging.root.level == logging.DEBUG:
                statusStr = '    Finished processing %d of %d (%.2f%%) bins.' % (
                    processed, total_bins, float(processed) * 100 / total_bins)
                sys.stdout.write('\r{}'.format(statusStr))
                sys.stdout.flush()

    def __run_prodigal(self, ttable):

        self.threads_per_bin = max(1, int(self.total_threads / len(self.bin_files)))
        logging.info("Calling genes in {} bins with {} threads:".format(len(self.bin_files), self.total_threads))

        # process each bin in parallel
        workerQueue = mp.Queue()
        writerQueue = mp.Queue()

        for bin in self.bin_files:
            workerQueue.put(bin)

        for _ in range(self.total_threads):
            workerQueue.put(None)

        used_ttables = mp.Manager().dict()
        coding_density = mp.Manager().dict()
        N50 = mp.Manager().dict()
        avg_gene_len = mp.Manager().dict()
        total_bases = mp.Manager().dict()
        cds_count = mp.Manager().dict()
        GC = mp.Manager().dict()
        totalContigs = mp.Manager().dict()
        maxContigLen = mp.Manager().dict()


        try:
            calcProc = []
            for _ in range(self.total_threads):
                calcProc.append(
                    mp.Process(target=self.__set_up_prodigal_thread, args=(workerQueue, writerQueue, ttable,
                                                                           used_ttables, coding_density,
                                                                           N50, avg_gene_len,
                                                                           total_bases, cds_count, GC, totalContigs, maxContigLen)))
            writeProc = mp.Process(target=self.__reportProgress, args=(len(self.bin_files), writerQueue))

            writeProc.start()

            for p in calcProc:
                p.start()

            for p in calcProc:
                p.join()

            writerQueue.put((None, None, None, None, None, None, None, None, None, None))
            writeProc.join()
        except:
            # make sure all processes are terminated
            for p in calcProc:
                p.terminate()

            writeProc.terminate()

        return used_ttables, coding_density, N50, avg_gene_len, total_bases, cds_count, GC, totalContigs, maxContigLen

    def __calculate_metadata(self, faa_files):

        self.threads_per_bin = max(1, int(self.total_threads / len(faa_files)))
        logging.info("Calculating metadata for {} bins with {} threads:".format(len(faa_files), self.total_threads))

        # process each bin in parallel
        workerQueue = mp.Queue()
        writerQueue = mp.Queue()

        for faa in faa_files:
            workerQueue.put(faa)

        for _ in range(self.total_threads):
            workerQueue.put(None)

        metadata_dict = mp.Manager().dict()

        try:
            calcProc = []
            for _ in range(self.total_threads):
                calcProc.append(
                    mp.Process(target=self.__set_up_metadata_thread, args=(workerQueue, writerQueue, metadata_dict)))
            writeProc = mp.Process(target=self.__report_progress_metadata, args=(len(faa_files), writerQueue))

            writeProc.start()

            for p in calcProc:
                p.start()

            for p in calcProc:
                p.join()

            writerQueue.put((None, None))
            writeProc.join()
        except:
            # make sure all processes are terminated
            for p in calcProc:
                p.terminate()

            writeProc.terminate()

        # metadata_dict = process into df (metadata_dict)
        return metadata_dict

    def __set_up_metadata_thread(self, queue_in, queue_out, metadata_dict):

        while True:
            bin = queue_in.get(block=True, timeout=None)
            if bin == None:
                break

            metadata_thread = metadata.MetadataCalculator(bin)
            name1, cdscount_series = metadata_thread.calculate_CDS()
            name2, aalength_series = metadata_thread.calculate_amino_acid_length()
            name3, aa_list, aa_counts = metadata_thread.calculate_amino_acid_counts()

            if name1 == name2 == name3:
                meta_thread_df = pd.DataFrame(
                    {'Name': [name1], 'CDS': [cdscount_series], 'AALength': [aalength_series]})
                for idx, aa in enumerate(aa_list):
                    meta_thread_df[aa] = aa_counts[idx]
            else:
                logging.error('Inconsistent name information in metadata calculation. Exiting.')
                sys.exit(1)

            metadata_dict[bin] = meta_thread_df

            queue_out.put(bin)

    def __report_progress_metadata(self, total_bins, queueIn):
        """Report number of processed bins."""

        processed = 0

        while True:
            bin = queueIn.get(block=True, timeout=None)
            if bin[0] == None:
                if logging.root.level == logging.INFO or logging.root.level == logging.DEBUG:
                    sys.stdout.write('\n')
                    sys.stdout.flush()
                break

            processed += 1

            if logging.root.level == logging.INFO or logging.root.level == logging.DEBUG:
                statusStr = '    Finished processing %d of %d (%.2f%%) bin metadata.' % (
                    processed, total_bins, float(processed) * 100 / total_bins)
                sys.stdout.write('\r{}'.format(statusStr))
                sys.stdout.flush()
