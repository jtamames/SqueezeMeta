#!/usr/bin/python3

"""
'-------------------------------------------------------------------------------'
Load results from SQM in anvi’o.
Run this script in a directory that contains:
    -Directory with results of sqm2anvio.pl.This directory will not be modified.
    It contains data prepared to be loaded into anvi'o:
        -Text files: contigs, genes, functions, taxonomy of contigs and genes and bins (optional).
        -Bam: directory with bam files. Bam files' names must finish with '-RAW.bam'.
    
USAGE: python3 anvi-load-sqm.py <-p project>
        [--num-threads int] [--run-HMMS] [--min-contigs-length int] [--min-mean-coverage int] [--skip-SNV-profiling] [--profile-SCVs]
    
IMPORTANT:
    This script do not overwrte anything.
    If you want to run it a second time, please remove all the previous created files.
    
    WARNING:
    Anvi’o is an advanced analysis and visualization platform for ‘omics data.
    Citation: Eren AM, Esen ÖC, Quince C, Vineis JH, Morrison HG, Sogin ML, Delmont TO. (2015)
              Anvi’o: an advanced analysis and visualization platform for ‘omics data. PeerJ 3:e1319
    This script is prepare to load data from SQM into anvi'o in a very straight way.
    If you want to explore more options of anvi’o please see the anvi'o project page
    for information, tutorials and more details. (http://merenlab.org/software/anvio/)
'-------------------------------------------------------------------------------'
"""

import subprocess
import os
import argparse
import sqlite3
from collections import defaultdict
from glob import glob
# Describe functions

#FUN: Parser input data
def parse_arguments():
    """Parse the command line arguments and return an object containing them."""
    # Required
    general = argparse.ArgumentParser(description='Process input files')
    general.add_argument('-p', '--project',type=str, required=True, help='Name of the directory product of sqm2anvio.pl')
    # Optional:
    general.add_argument('--num-threads', default=12, type=int, help='Number of threads.')
    general.add_argument('--run-HMMS', action = 'store_true', help='Run all 4 HMM profiles (anvi’o function). These are the currently available ones in anvi’o: "Rinke_et_al" , "Campbell_et_al" ,"BUSCO_83_Protista" , "Ribosomal_RNAs".')
    general.add_argument('--min-contigs-length', default=0, type=int, help='Anvi’o parameter: Minimum length of contigs in a BAM file to analyze. The minimum length should be long enough for tetranucleotide frequency analysis to be meaningful.')  
    general.add_argument('--min-mean-coverage', default=0, type=int, help='Anvi’o parameter: Minimum mean coverage for contigs to be kept in the analysis. The default value is 0, which is for your best interest if you are going to profile multiple BAM files which are then going to be merged for a cross sectional or time series analysis.')
    general.add_argument('--skip-SNV-profiling', action='store_true', help='The use of this flag will instruct profiler to skip that step.')
    general.add_argument('--profile-SCVs', action='store_true', help='Anvi’o parameter to perform accurate characterization of codon frequencies in genes during profiling.')
    
    args = general.parse_args()
    
    return(args)

#FUN: Run and check whether the subprocess worked well
def run_command(command):# run_command
    """Run the command and check the success of the subprocess. Return exit if it went wrong."""
    exitcode=subprocess.call(map(str,command)) 
    if exitcode!=0:
        print('There must be some problem with {}.\nIt\'s better to stop and check it.'.format(' '.join(command)))
        exit(-1)

#FUN: Check arguments
def check_argument(project):
    #Check if the directory exists and contains at least the required arguments:
    if os.path.isdir('./{}'.format(project)):
        files=[]
        for f in ['_anvio_contigs', '_anvio_genes', '_anvio_functions','_anvio_taxonomy','_anvio_contig_taxonomy']:
            globList = glob('{}/*{}.txt'.format(project,f)) #recover names of files that are arguments for anvio
            if not globList:
                print('The file ended in \'{}.txt\' seems not to exist. This file should have been generated by sqm2qnvio.pl'.format(f))
                exit(-1)
            elif len(globList) > 1:
                print('There seems to be more than one file ended in \'{}.txt\'. Are you trying to trick us?'.format(f))
                exit(-1)
            elif not os.path.isfile(globList[0]):
                print('{} is not a file! We really thought this would never happen!'.format(globlist[0]))
                exit(-1)
            else:
                files.append(globList[0])
        
        if not os.path.isdir('./{}/bam'.format(project)):
            print('There is not a bam directory. It should have been generated by sqm2anvio.pl')
            exit(-1)
        else: files.append('bam')           

        return files
    else:
        print('{} does not exist in'.format(project,os.getcwd()))
        exit(-1)

#FUN: Main
def main(args):
    """Get things done"""
    
    print('------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------')
    print('Loading results from SQM in anvi’o.\nAnvi’o is an advanced analysis and visualization platform for ‘omics data.',
          'Citation: Eren AM, Esen ÖC, Quince C, Vineis JH, Morrison HG, Sogin ML, Delmont TO. (2015) Anvi’o: an advanced analysis and visualization platform for ‘omics data. PeerJ 3:e1319',
          'This script is prepare to load data from SQM into anvi’o in a very straight way.',
          'If you want to explore more options of anvi’o',
          'please consider to see the anvi’o project page for information, tutorials and more details. (http://merenlab.org/software/anvio/)', sep='\n')
    print('------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------')
    

    #! FIRST: check if the anvi'o enviroment is open
    # If the process is open, it'll return 0, otherwise not
    exitcode = subprocess.call(['which', 'anvi-self-test'],stdout=subprocess.PIPE)
    if exitcode!=0:
        print('Anvi’o has not been detected. Are you sure that it has been activated?')
        exit(-1)
    else:
        contigs, genes, functions, tax_genes, tax_contigs, bams = check_argument(args.project)
        #! SECOND: Start running anvio, independently of number of samples. Make complete CONTIGS.db
        print('Anvi’o is preparing your contig database: loading contigs, genes, functions and taxonomy!')
        # Load contigs & genes with some parameters from anvio
        command=['anvi-gen-contigs-database', '-f', contigs,'-n',args.project,'-o', 'CONTIGS.db',  '--external-gene-call', genes,'--ignore-internal-stop-codons']
        run_command(command)
        print('Your contigs database is CONTIGS.db. Contigs and genes have been loaded')
        
        # Run if it's required the HMMS option from anvio
        if args.run_HMMS:
            print('You selected to run HMMS to detect SCGs')
            command=['anvi-run-hmms', '-c', 'CONTIGS.db', '--num-threads' ,args.num_threads]
            run_command(command)
        
        # Load functions    
        command=['anvi-import-functions', '-c', 'CONTIGS.db', '-i' , functions]
        run_command(command)
        print('Your functions have been loaded')
        
        #Load taxonomy
        conn = sqlite3.connect('CONTIGS.db')
        c = conn.cursor()
        
        ## Execute selection of splits
        c.execute('SELECT split FROM splits_basic_info')
        splits = [s[0] for s in c.fetchall()]
        
        ## Dictionary: {contig:tax }
        with open(tax_contigs,'r') as infile:
            contig_tax={}
            infile.readline() # burn headers
            contig_tax = dict([line.rstrip('\n').split('\t',1) for line in infile]) #remove superkingdom
        
        split_tax = {s: contig_tax[s.split('_split_')[0]] for s in splits}
        
        
        ## Dictionary {gen:tax}
        with open(tax_genes,'r') as infile2:
            gen_tax={}
            infile2.readline() # burn headers
            gen_tax = dict([line.rstrip('\n').split('\t',1) for line in infile2])
                
        
        ## Join all the possible taxonomy combinations: from genes + contigs # could happen that genes have different taxonomy annotations than in contigs are not consiodered due to the consensus procedure.
        taxa = set( list(split_tax.values()) + list(gen_tax.values()) )
        taxon_names = {taxon:i+1 for i,taxon in enumerate(taxa)}
        
        
        for t in taxon_names.keys():
            
            stmt = "INSERT INTO taxon_names (taxon_id, t_phylum, t_class, t_order, t_family, t_genus, t_species) VALUES (?,?,?,?,?,?,?);"
            tax_list=t.split('\t')[1:len(t)]
            values=[taxon_names[t]]+tax_list #remove superkingdom
            c.execute(stmt,values)
            
        #splits_taxonomy: split \t ids de la taxonomy   
        for s in split_tax:
            stmt="INSERT INTO splits_taxonomy (split, taxon_id) VALUES (?,?);"
            values=[s,taxon_names[split_tax[s]]]
            c.execute(stmt,values)
            
        #genes_taxonomy: gen \t ids de la taxonomy   
        for gen in gen_tax:
            stmt="INSERT INTO genes_taxonomy (gene_callers_id, taxon_id) VALUES (?,?);"
            values=[gen,taxon_names[gen_tax[gen]]]
            c.execute(stmt,values)
        #Close DB saving changes    
        conn.commit()
        print('Your taxonomy has been loaded')
                        
            
    #! THIRD: Make profiles.db. Load bam files, important: distinguish number of samples
        print('Anvi\'o will load your bam files')
        samples = [f for f in os.listdir('{}/bam'.format(args.project)) if f.endswith('-RAW.bam')]
        profiles_names= ['{}_{}_{}'.format(f.replace('-RAW.bam',''),args.min_contigs_length,args.min_mean_coverage) for f in samples]
        if not samples:
                profile_name='BLANK_PROFILE_{}'.format(args.mincontigs_length)
                command_base=['anvi-profile', '-o', profile_name, '-c', 'CONTIGS.db', '--blank-profile','-S','Blank','--min-contig-length',args.min_contigs_length ,'--num-threads' ,args.num_threads,'--skip-hierarchical-clustering']
                run_command(command_base)
                #Load Bin collection
                if len(glob('{}/*anvio_bins.txt'.format(args.project)))==1:
                    if os.path.isfile(glob('{}/*anvio_bins.txt'.format(args.project))[0]):
                        print('Anvi’o will load your DAS collection')
                        samples_profiles=['{}/PROFILE.db'.format(profile_name)]
                        command=['anvi-import-collection',glob('{}/*anvio_bins.txt'.format(args.project))[0],'-c','CONTIGS.db','-p'] + samples_profiles + ['-C', 'DAS', '--contigs-mode']
                        run_command(command)
                        print('Anvi’o finished to load your DAS collection')
                    else: print('Your bins file is not a file! We really thought this would never happen!')
                else:print('There is none bins collection')
                
        else:
            # Iterate over directory with bam files
            for f in samples:
                #print(f)
                print('Processing {}'.format(f))
                f_in=args.project + '/bam/' + f
                f_out=f_in.replace('-RAW','')
                # Order & Index bam file
                command=['anvi-init-bam', f_in, '-o', f_out]
                run_command(command)
                #print(f_out)
                print('Profiling {}'.format(f))
                # Profile bam file
                
                profile_name='{}_{}_{}'.format(f.replace('-RAW.bam',''),args.min_contigs_length,args.min_mean_coverage)
                command_base=['anvi-profile','-i',f_out, '-o', profile_name, '-c', 'CONTIGS.db', '--min-contig-length', args.min_contigs_length ,'--min-mean-coverage', args.min_mean_coverage,'--num-threads' ,args.num_threads,'--skip-hierarchical-clustering']
                if args.skip_SNV_profiling:
                    command_base.append('--skip-SNV-profiling')
                if args.profile_SCVs:
                    command_base.append('--profile-SCVs')
                if len(samples)==1:
                    run_command(command_base)
                    print('Removing sort and index bam file')
                    command=['rm',f_out,'{}.bai'.format(f_out) ]
                    run_command(command)    
                    
                    #Load Bin collection
                    if len(glob('{}/*anvio_bins.txt'.format(args.project)))==1:
                        if os.path.isfile(glob('{}/*anvio_bins.txt'.format(args.project))[0]):
                            print('Anvi\'o will load your DAS collection')
                            samples_profiles=['{}/PROFILE.db'.format(profile_name)]
                            command=['anvi-import-collection',glob('{}/*anvio_bins.txt'.format(args.project))[0],'-c','CONTIGS.db','-p'] + samples_profiles + ['-C', 'DAS', '--contigs-mode']
                            run_command(command)
                            print('Anvi’o finished to load your DAS collection')
                        else: print('Your bins file is not a file! We really thought this would never happen!')
                    else:print('There is none bins collection')
                    
                    
                else:
                    run_command(command_base)
                    print('Removing sort and index bam file')
                    command=['rm',f_out, '{}.bai'.format(f_out)]
                    run_command(command)    
                    
                
    # If there are 2 or more samples: Merge profiles.db
        if len(samples) > 1:
            print('Anvi’o will merge your profile databases')
            profile_dir=['{}/PROFILE.db'.format(p) for p in profiles_names]
            print('These are the profiles db {} that will be merged'.format(profile_dir))
            command=['anvi-merge'] + profile_dir + ['-c','CONTIGS.db','-o','SAMPLES-MERGED_{}_{}'.format(args.min_contigs_length,args.min_mean_coverage),'--skip-concoct-binning','--skip-hierarchical-clustering']
            run_command(command)
        
            #Load Bin collection
            if len(glob('{}/*anvio_bins.txt'.format(args.project)))==1:
                if os.path.isfile(glob('{}/*anvio_bins.txt'.format(args.project))[0]):
                    print('Anvi’o will load your DAS collection')
                    samples_profiles=['{}/PROFILE.db'.format(profile_name)]
                    command=['anvi-import-collection',glob('{}/*anvio_bins.txt'.format(args.project))[0],'-c','CONTIGS.db','-p', 'SAMPLES-MERGED_{}_{}/PROFILE.db'.format(args.min_contigs_length,args.min_mean_coverage), '-C', 'DAS', '--contigs-mode']
                    run_command(command)
                    print('Anvi’o finished to load your DAS collection')
                else: print('Your bins file is not a file! We really thought this would never happen!')
            else:print('There is no bins collection. Skipping...')
        else:
            print('Anvi’o finished to load your data :)')
        
################################################################################################################
    
if __name__ == '__main__':
    main(parse_arguments())