#!/usr/bin/python

import sqlite3
from subprocess import call



def subset_anvio(splits, contigs_db, profile_db, outdir):

    splits = list(splits)
    splitsSub = '({})'.format(','.join('?' * len(splits)))
    contigs = [split.rsplit('_split', 1)[0] for split in splits]
    contigsSub = '({})'.format(','.join('?' * len(contigs)))
    
    conn = sqlite3.connect(contigs_db)
    c = conn.cursor()
    query = 'SELECT gene_callers_id FROM genes_in_splits WHERE split IN {}'.format(splitsSub)
    c.execute(query, splits)
    genes =  [x[0] for x in c.fetchall()]
    conn.close()
    genesSub = '({})'.format(','.join('?' * len(genes)))

    aux_db = '{}/AUXILIARY-DATA.db'.format('/'.join(profile_db.split('/')[:-1]))
    new_contigs_db = '{}/CONTIGS.db'.format(outdir)
    new_profile_db = '{}/PROFILE.db'.format(outdir)
    new_aux_db = '{}/AUXILIARY-DATA.db'.format(outdir)
    
    call(['mkdir', outdir])

    max_variable_limit = max([len(splits), len(genes)])

    ### Filter contigs.db
    print('Filtering {}'.format(contigs_db))
    conn = sqlite3.connect(new_contigs_db)
    c = conn.cursor()
    c.execute('ATTACH DATABASE ? AS old', (contigs_db,) )

    # self
    c.execute('CREATE TABLE self AS SELECT * FROM old.self')
    c.execute('UPDATE self SET value=? WHERE key="num_contigs";', (len(contigs),) )
    c.execute('UPDATE self SET value=? WHERE key="num_splits";', (len(splits),) )
    #UPDATE self SET value=2 WHERE key='total_length'; # calculate total length!
    c.execute('UPDATE self SET value="fakehash" WHERE key="contigs_db_hash";')
    # hmm_hits
    c.execute('CREATE TABLE hmm_hits AS SELECT * FROM old.hmm_hits') # Just copy
    # hmm_hits_info
    c.execute('CREATE TABLE hmm_hits_info AS SELECT * FROM old.hmm_hits_info') # Just copy
    # hmm_hits_in_splits
    query = 'CREATE TABLE hmm_hits_in_splits AS SELECT * FROM old.hmm_hits_in_splits WHERE split IN {}'.format(splitsSub)
    c.execute(query, splits)
    # collections_info
    c.execute('CREATE TABLE collections_info AS SELECT * FROM old.collections_info') # Just copy
    # collections_bins_info
    c.execute('CREATE TABLE collections_bins_info AS SELECT * FROM old.collections_bins_info') # Just copy
    # collections_of_contigs
    query = 'CREATE TABLE collections_of_contigs AS SELECT * FROM old.collections_of_contigs WHERE contig IN {}'.format(contigsSub)
    c.execute(query, contigs)
    # collections_of_splits
    query = 'CREATE TABLE collections_of_splits AS SELECT * FROM old.collections_of_splits WHERE split IN {}'.format(splitsSub)
    c.execute(query, splits)
    # genes_in_contigs
    query = 'CREATE TABLE genes_in_contigs AS SELECT * FROM old.genes_in_contigs WHERE contig IN {};'.format(contigsSub)
    c.execute(query, contigs)
    # genes_in_splits
    query = 'CREATE TABLE genes_in_splits AS SELECT * FROM old.genes_in_splits WHERE split IN {};'.format(splitsSub)
    c.execute(query, splits)
    # splits_taxonomy
    query = 'CREATE TABLE splits_taxonomy AS SELECT * FROM old.splits_taxonomy WHERE split IN {}'.format(splitsSub)
    c.execute(query, splits)
    # taxon_names
    c.execute('CREATE TABLE taxon_names AS SELECT * FROM old.taxon_names') # Just copy
    # genes_taxonomy
    query = 'CREATE TABLE genes_taxonomy AS SELECT * FROM old.genes_taxonomy WHERE gene_callers_id IN {}'.format(genesSub)
    c.execute(query, genes)
    # contig_sequences
    query = 'CREATE TABLE contig_sequences AS SELECT * FROM old.contig_sequences WHERE contig IN {}'.format(contigsSub)
    c.execute(query, contigs)
    # gene_functions
    query = 'CREATE TABLE gene_functions AS SELECT * FROM old.gene_functions WHERE gene_callers_id IN {}'.format(genesSub)
    c.execute(query, genes)
    # gene_amino_acid_sequences
    query = 'CREATE TABLE gene_amino_acid_sequences AS SELECT * FROM old.gene_amino_acid_sequences WHERE gene_callers_id IN {}'.format(genesSub)
    c.execute(query, genes)
    # splits_basic_info
    query = 'CREATE TABLE splits_basic_info AS SELECT * FROM old.splits_basic_info WHERE split IN {}'.format(splitsSub)
    c.execute(query, splits)
    # contigs_basic_info
    query = 'CREATE TABLE contigs_basic_info AS SELECT * FROM old.contigs_basic_info WHERE contig IN {}'.format(contigsSub)
    c.execute(query, contigs)
    # nt_position_info
    query = 'CREATE TABLE nt_position_info AS SELECT * FROM old.nt_position_info WHERE contig_name IN {}'.format(contigsSub)
    c.execute(query, contigs)
    # kmer_contigs
    query = 'CREATE TABLE kmer_contigs AS SELECT * FROM old.kmer_contigs WHERE contig IN {}'.format(splitsSub) # table is called kmer_contigs but they seem to be using splits as indexes.
    c.execute(query, splits)
    # kmer_splits
    query = 'CREATE TABLE kmer_splits AS SELECT * FROM old.kmer_splits WHERE contig IN {}'.format(splitsSub)
    c.execute(query, splits)
    
    conn.commit()
    conn.close()
    #############################################################################################

    ### Filter profile.db
    print('Filtering {}'.format(profile_db))
    conn = sqlite3.connect(new_profile_db)
    c = conn.cursor()
    c.execute('ATTACH DATABASE ? AS old', (profile_db,) )

    # self
    c.execute('CREATE TABLE self AS SELECT * FROM old.self')
    c.execute('UPDATE self SET value=? WHERE key="num_contigs";', (len(contigs),) )
    c.execute('UPDATE self SET value=? WHERE key="num_splits";', (len(splits),) )
    # UPDATE TOTAL LENGTH AND TOTAL READS MAPPED AND WHATNOT
    c.execute('UPDATE self SET value="fakehash" WHERE key="contigs_db_hash";')
    # item_additional_data
    c.execute('CREATE TABLE item_additional_data AS SELECT * FROM old.item_additional_data WHERE NULL') # empty
    # item_orders
    c.execute('CREATE TABLE item_orders AS SELECT * FROM old.item_orders WHERE NULL') # empty
    # layer_additional_data
    c.execute('CREATE TABLE layer_additional_data AS SELECT * FROM old.layer_additional_data WHERE NULL') # empty
    # layer_orders
    c.execute('CREATE TABLE layer_orders AS SELECT * FROM old.layer_orders WHERE NULL') # empty
    # variable_nucleotides
    query = 'CREATE TABLE variable_nucleotides AS SELECT * FROM old.variable_nucleotides WHERE split_name NOT IN {}'.format(splitsSub)
    c.execute(query, splits)
    # variable_codons
    query = 'CREATE TABLE variable_codons AS SELECT * FROM old.variable_codons WHERE corresponding_gene_call IN {}'.format(genesSub)
    c.execute (query, genes)
    # views
    c.execute('CREATE TABLE views AS SELECT * FROM old.views') # Just copy
    # collections_info
    c.execute('CREATE TABLE collections_info AS SELECT * FROM old.collections_info') # Just copy
    # collections_bins_info
    c.execute('CREATE TABLE collections_bins_info AS SELECT * FROM old.collections_bins_info') # Just copy
    # collections_of_contigs
    c.execute('CREATE TABLE collections_of_contigs AS SELECT * FROM old.collections_of_contigs') # Just copy
    # collections_of_splits
    c.execute('CREATE TABLE collections_of_splits AS SELECT * FROM old.collections_of_splits') # Just copy
    # states
    c.execute('CREATE TABLE states AS SELECT * FROM old.states') # Just copy
    # std_coverage_contigs
    query = 'CREATE TABLE std_coverage_contigs AS SELECT * FROM old.std_coverage_contigs WHERE contig IN {}'.format(splitsSub)
    c.execute(query, splits)
    # mean_coverage_contigs
    query = 'CREATE TABLE mean_coverage_contigs AS SELECT * FROM old.mean_coverage_contigs WHERE contig IN {}'.format(splitsSub)
    c.execute(query, splits)
    # mean_coverage_Q2Q3_contigs
    query = 'CREATE TABLE mean_coverage_Q2Q3_contigs AS SELECT * FROM old.mean_coverage_Q2Q3_contigs WHERE contig IN {}'.format(splitsSub)
    c.execute(query, splits)
    # max_normalized_ratio_contigs
    query = 'CREATE TABLE max_normalized_ratio_contigs AS SELECT * FROM old.max_normalized_ratio_contigs WHERE contig IN {}'.format(splitsSub)
    c.execute(query, splits)
    # relative_abundance_contigs
    query = 'CREATE TABLE relative_abundance_contigs AS SELECT * FROM old.relative_abundance_contigs WHERE contig IN {}'.format(splitsSub)
    c.execute(query, splits)
    # detection_contigs
    query = 'CREATE TABLE detection_contigs AS SELECT * FROM old.detection_contigs WHERE contig IN {}'.format(splitsSub)
    c.execute(query, splits)    
    # abundance_contigs
    query = 'CREATE TABLE abundance_contigs AS SELECT * FROM old.abundance_contigs WHERE contig IN {}'.format(splitsSub)
    c.execute(query, splits)    
    # variability_contigs
    query = 'CREATE TABLE variability_contigs AS SELECT * FROM old.variability_contigs WHERE contig IN {}'.format(splitsSub)
    c.execute(query, splits)    
    # std_coverage_splits
    query = 'CREATE TABLE std_coverage_splits AS SELECT * FROM old.std_coverage_splits WHERE contig IN {}'.format(splitsSub)
    c.execute(query, splits)
    # mean_coverage_splits
    query = 'CREATE TABLE mean_coverage_splits AS SELECT * FROM old.mean_coverage_splits WHERE contig IN {}'.format(splitsSub)
    c.execute(query, splits)
    # mean_coverage_Q2Q3_splits
    query = 'CREATE TABLE mean_coverage_Q2Q3_splits AS SELECT * FROM old.mean_coverage_Q2Q3_splits WHERE contig IN {}'.format(splitsSub)
    c.execute(query, splits)    
    # max_normalized_ratio_splits
    query = 'CREATE TABLE max_normalized_ratio_splits AS SELECT * FROM old.max_normalized_ratio_splits WHERE contig IN {}'.format(splitsSub)
    c.execute(query, splits)    
    # relative_abundance_splits
    query = 'CREATE TABLE relative_abundance_splits AS SELECT * FROM old.relative_abundance_splits WHERE contig IN {}'.format(splitsSub)
    c.execute(query, splits)
    # detection_splits
    query = 'CREATE TABLE detection_splits AS SELECT * FROM old.detection_splits WHERE contig IN {}'.format(splitsSub)
    c.execute(query, splits)
    # abundance_splits
    query = 'CREATE TABLE abundance_splits AS SELECT * FROM old.abundance_splits WHERE contig IN {}'.format(splitsSub)
    c.execute(query, splits)
    # variability_splits
    query = 'CREATE TABLE variability_splits AS SELECT * FROM old.variability_splits WHERE contig IN {}'.format(splitsSub)
    c.execute(query, splits)

    conn.commit()
    conn.close()

    #############################################################################################

    ### Filter aux.db
    print('Filtering {}'.format(aux_db))
    conn = sqlite3.connect(new_aux_db)
    c = conn.cursor()
    c.execute('ATTACH DATABASE ? AS old', (aux_db,) )

    # self
    c.execute('CREATE TABLE self AS SELECT * FROM old.self')
    c.execute('UPDATE self SET value="fakehash" WHERE key="contigs_db_hash";')
    # split_coverages
    query = 'CREATE TABLE split_coverages AS SELECT * FROM old.split_coverages WHERE split_name IN {}'.format(splitsSub)
    c.execute(query, splits)

    conn.commit()
    conn.close()



#subset_anvio(['Merged_10000_split_00001', 'Merged_47_split_00001', 'Merged_100_split_00001'], '/media/disk8/fer/Hadza/contigs.db', '/media/disk8/fer/Hadza/SAMPLES-MERGED/PROFILE.db', 'test')
