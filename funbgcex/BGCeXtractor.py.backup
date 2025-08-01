#!/usr/bin/env python3

import logging
import os
from pathlib import Path
import pandas as pd
from natsort import natsorted
import shutil
import glob
import sys
import time
import multiprocessing
from multiprocessing import Pool, Manager, cpu_count, Lock
import psutil
from tqdm import tqdm
from funbgcex.GeneralCommands import *
from funbgcex.SeqHandling import *
from funbgcex.BGfinder import *
from funbgcex.HTMLgenerator import HTMLgenerator
from funbgcex.SimilarBGCfinder import MakeProtBGCidDict
from funbgcex.MemoryOptimizer import MemoryMonitor


# Using the enhanced MemoryMonitor from MemoryOptimizer module


def get_optimal_worker_count(requested_workers=None, max_memory_gb=None, min_workers=1, reserve_cores=1):
    """
    Determine the optimal number of worker processes based on system resources.
    
    Args:
        requested_workers: User-specified number of workers (None for auto-detect)
        max_memory_gb: Maximum memory usage in GB (None for no limit)
        min_workers: Minimum number of workers to use
        reserve_cores: Number of CPU cores to reserve for system processes
    
    Returns:
        Optimal number of worker processes
    """
    # Get available CPU cores
    available_cores = cpu_count()
    
    # Start with requested workers or auto-detect
    if requested_workers is None:
        # Auto-detect: leave at least reserve_cores for system
        workers = max(min_workers, available_cores - reserve_cores)
    else:
        # Use user-specified value but ensure it's within valid range
        workers = max(min_workers, min(requested_workers, available_cores))
    
    # Adjust for memory constraints if specified
    if max_memory_gb is not None:
        try:
            # Get available memory in GB
            available_memory = psutil.virtual_memory().available / (1024**3)
            
            # Estimate memory usage per worker (conservative estimate: 2GB per worker)
            memory_per_worker = 2.0
            
            # Calculate maximum workers based on memory
            max_workers_by_memory = int(available_memory / memory_per_worker)
            
            # Use the minimum of CPU-based and memory-based worker counts
            workers = min(workers, max_workers_by_memory)
            
            if workers < min_workers:
                print(f"Warning: Limited memory reduces worker count from {requested_workers or available_cores - reserve_cores} to {workers}")
        except Exception as e:
            print(f"Warning: Could not determine memory constraints: {e}")
    
    return workers


def retry_operation(operation, args, max_retries=3, delay=1):
    """
    Retry an operation with exponential backoff if it fails.
    
    Args:
        operation: Function to retry
        args: Arguments to pass to the function
        max_retries: Maximum number of retry attempts
        delay: Initial delay between retries in seconds
    
    Returns:
        Result of the operation or None if all retries failed
    """
    for attempt in range(max_retries):
        try:
            return operation(*args)
        except Exception as e:
            if attempt == max_retries - 1:
                # Last attempt failed, return None
                return None
            # Exponential backoff
            wait_time = delay * (2 ** attempt)
            time.sleep(wait_time)
    return None


def process_single_file_worker(args):
    """
    Worker function for processing a single GenBank file in parallel.
    
    Args:
        args: Tuple containing (file, gbk_dir, results_dir, mode, query, gap_allowed, 
              max_bgc_gap, min_prot_len, num_of_genes_checked, min_identity, 
              current_dir, IDdict, GeneNumDict, MetabDict, all_cluster_csv)
    
    Returns:
        Tuple containing (file_name, BGC_count, output_dir, success, error_message)
    """
    file, gbk_dir, results_dir, mode, query, gap_allowed, max_bgc_gap, min_prot_len, num_of_genes_checked, min_identity, current_dir, IDdict, GeneNumDict, MetabDict, all_cluster_csv, csv_lock, directory_lock = args
    
    try:
        file_name = Path(file).stem
        output_dir = f"{results_dir}/{file_name}_results"
        # Create unique temporary directory with process ID to avoid conflicts
        import os
        process_id = os.getpid()
        temp_dir = f"{output_dir}/temp_{process_id}"
        temp_dir_ = f"{temp_dir}/temp"
        os.makedirs(output_dir, exist_ok=True)
        os.makedirs(temp_dir, exist_ok=True)
        df = pd.DataFrame()
        
        # Create a file-specific logger
        worker_logger = logging.getLogger(f"worker_{file_name}")
        worker_logger.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s - %(levelname)s:%(name)s - %(message)s')
        log_file = logging.FileHandler(f'{output_dir}/{file_name}.log')
        log_file.setFormatter(formatter)
        worker_logger.addHandler(log_file)
        
        worker_logger.info(f"{file_name}: BGC extraction started")
        
        # Make a file to which the information of all clusters will be saved
        columns = ["BGC no.","Scaffold","Start position","End position","Locus tag start","Locus tag end",
                    "Number of genes","Core enzymes","Similar BGC","Similarity score","Metabolite from similar BGC","Pfam domains"]
        
        worker_logger.info(f"{mode} mode was selected")
        
        BGCdf = pd.DataFrame(columns=columns)
        cluster_csv = f"{output_dir}/BGCs.csv"
        BGCdf.to_csv(cluster_csv)
        
        # Check the input file. If locus tags are not given for CDSs, they are added as PROTEINXXXXXX.
        gbk_file = SeqHandling(file, df)
        gbk_file.AddLocusTag(temp_dir)
        modified_gbk_file = glob.glob(f"{temp_dir}/*.gbk")[0]
        modified_gbk = SeqHandling(modified_gbk_file, df)
        
        worker_logger.debug("Processed the input GenBank file")
        
        # Make necessary directories
        CDS_dir = f"{temp_dir}/extracted_CDS"
        GBK_dir = f"{temp_dir}/scaffolds"
        
        os.makedirs(CDS_dir, exist_ok=True)
        os.makedirs(GBK_dir, exist_ok=True)
        
        # Obtain the fungus's name
        fungus_name = modified_gbk.name()  # Special symbols that might cause errors removed
        fungus_name_original = modified_gbk.name_o()  # Original fungus's name
        
        worker_logger.debug(f"Original fungus's name: {fungus_name_original}")
        worker_logger.debug(f"Modified version of fungus's name: {fungus_name}")
        
        # Extract all CDSs as a single fasta file and each scaffold as a gbk file
        modified_gbk.CDS_extractor(CDS_dir)
        modified_gbk.GBK_extractor(GBK_dir)
        fasta_file = f"{CDS_dir}/{file_name}.fasta"
        
        worker_logger.debug("All protein sequences were extracted")
        worker_logger.debug("All scaffolds were extracted")
        
        BGC_count = 0
        
        if os.stat(fasta_file).st_size != 0:
            # Create unique BLAST database directory for each worker process
            dmnd_db_dir = f"{temp_dir}/BLASTdb_{process_id}"
            os.makedirs(dmnd_db_dir, exist_ok=True)
            
            if mode == "all":
                # hmmscan to detect core enzymes
                hmmscan_output_dir = f"{temp_dir}/hmmscan_{process_id}"
                os.makedirs(hmmscan_output_dir, exist_ok=True)
                hmmscan_result = f"{hmmscan_output_dir}/hmm_core.txt"
                # Create a copy of the database for each worker to avoid file locking issues
                import shutil
                hmmscan_core_database_copy = f"{temp_dir}/core_UstY_{process_id}.hmm"
                shutil.copy2(f"{current_dir}/data/hmm/core/core_UstY.hmm", hmmscan_core_database_copy)
                
                worker_logger.debug("Running hmmscan to detect core enzymes")
                # Use parallel hmmscan for better performance
                runHMMscan_parallel([fasta_file], hmmscan_core_database_copy, hmmscan_output_dir, "1e-5", workers=2)
                
                worker_logger.debug("Finding core enzymes")
                Core_extractor(hmmscan_result, df, temp_dir_)
            
            if mode == "target":
                # Create DIAMOND database with unique name
                target_dmnd_db = f"{dmnd_db_dir}/target_{process_id}"
                worker_logger.debug("Creating DIAMOND database using the input fasta file")
                makeDIAMONDdb(query, target_dmnd_db)
                
                # Perform blastp to find target homologues
                blastp_result_dir = f"{temp_dir}/blastp_results_{process_id}"
                os.makedirs(blastp_result_dir, exist_ok=True)
                
                target_blastp_result = f"{blastp_result_dir}/target_blastp_result.xml"
                worker_logger.debug("Running DIAMOND to find target proteins")
                # Use parallel DIAMOND for better performance
                RunDIAMOND_parallel([fasta_file], target_dmnd_db, blastp_result_dir, 1, 1, workers=2)
                
                # Add target homologue information to the csv file
                worker_logger.debug("Adding the information of the closest homologue of each target")
                AddTargetHomologue(target_blastp_result, df)
            
            if mode == "pfam":
                # Extract the target Pfam HMM file
                funSMhmm = f"{current_dir}/data/hmm/fungalSM/fungalSM.hmm"
                hmm_dir = f"{temp_dir}/hmm_{process_id}"
                os.makedirs(hmm_dir, exist_ok=True)
                target_hmm = f"{hmm_dir}/{query}.hmm"
                
                worker_logger.debug("Fetching the required HMM file")
                runHMMfetch(funSMhmm, query, target_hmm)
                
                if os.stat(target_hmm).st_size == 0:
                    worker_logger.error("The given protein family name was not found in the database.")
                    return (file_name, 0, output_dir, False, "The given protein family name was not found in the database.")
                
                worker_logger.debug("Running hmmpress")
                runHMMpress(target_hmm)
                
                # hmmscan to detect target proteins
                hmmscan_output_dir = f"{temp_dir}/hmmscan_pfam_{process_id}"
                os.makedirs(hmmscan_output_dir, exist_ok=True)
                hmmscan_result = f"{hmmscan_output_dir}/hmm_target.txt"
                hmmscan_target_database = f"{temp_dir}/hmm/{query}.hmm"
                
                worker_logger.debug("Running hmmscan to detect target proteins")
                # Use parallel hmmscan for better performance
                runHMMscan_parallel([fasta_file], hmmscan_target_database, hmmscan_output_dir, "1e-5", workers=2)
                
                worker_logger.debug("Adding target information")
                AddTargetPfam(hmmscan_result, df)
            
            # Extract protein sequences around core/target proteins
            extract_CDS_target_dir = f"{temp_dir}/extracted_CDS_target_{process_id}"
            os.makedirs(extract_CDS_target_dir, exist_ok=True)
            
            worker_logger.debug("Extracting protein sequences around core/target proteins")
            ExtractCDS4Check(mode, num_of_genes_checked, extract_CDS_target_dir, df)
            
            # Combine all extracted fasta files as a single file
            fasta_files = glob.glob(f"{extract_CDS_target_dir}/*.fasta")
            if len(fasta_files) != 0:
                worker_logger.debug("Combining all fasta files")
                combineFASTA(extract_CDS_target_dir, temp_dir)
                combined_fasta = f"{temp_dir}/combined.fasta"
                
                if mode == "target" or mode == "pfam":
                    # hmmscan to detect core enzymes
                    hmmscan_output_dir = f"{temp_dir}/hmmscan_core_{process_id}"
                    os.makedirs(hmmscan_output_dir, exist_ok=True)
                    hmmscan_result = f"{hmmscan_output_dir}/hmm_core.txt"
                    # Create a copy of the database for each worker to avoid file locking issues
                    hmmscan_core_database_copy = f"{temp_dir}/core_{process_id}.hmm"
                    shutil.copy2(f"{current_dir}/data/hmm/core/core.hmm", hmmscan_core_database_copy)
                    
                    worker_logger.debug("Running hmmscan to detect core enzymes")
                    # Use parallel hmmscan for better performance
                    runHMMscan_parallel([combined_fasta], hmmscan_core_database_copy, hmmscan_output_dir, "1e-5", workers=2)
                    
                    worker_logger.debug("Finding core enzymes")
                    Core_extractor(hmmscan_result, df, temp_dir_)
                
                worker_logger.debug("Finding potential RiPP precursor proteins")
                aa_len, unique_aa, min_repeat, min_match = 10, 3, 3, 9
                RiPPppFinder(aa_len, unique_aa, min_repeat, min_match, combined_fasta, df)
                
                # Perform hmmscan using the fungal SM protein domain database
                hmmscan_result_ = f"{hmmscan_output_dir}/hmmscan_result.txt"
                # Create a copy of the database for each worker to avoid file locking issues
                hmmscan_SM_database_copy = f"{temp_dir}/fungalSM_{process_id}.hmm"
                shutil.copy2(f"{current_dir}/data/hmm/fungalSM/fungalSM.hmm", hmmscan_SM_database_copy)
                
                worker_logger.debug("Running hmmscan to detect potential fungal SM proteins")
                # Use parallel hmmscan for better performance
                runHMMscan_parallel([combined_fasta], hmmscan_SM_database_copy, hmmscan_output_dir, "1e-5", workers=2)
                
                # Add Pfam domain information to the csv file
                SMhmm = f"{current_dir}/data/hmm/fungalSM/fungalSM.hmm"
                
                worker_logger.debug("Adding Pfam information")
                AddPfam(hmmscan_result_, df, SMhmm)
                
                # Find the closest homologue in the fungal biosynthetic proteins
                # Create a copy of the database for each worker to avoid file locking issues
                dmnd_db_SM_copy = f"{temp_dir}/fungal_SM_proteins_{process_id}"
                shutil.copytree(f"{current_dir}/data/DIAMOND/fungal_SM_proteins.dmnd", f"{dmnd_db_SM_copy}.dmnd")
                
                blastp_result_dir = f"{temp_dir}/blastp_results_SM_{process_id}"
                os.makedirs(blastp_result_dir, exist_ok=True)
                blastp_result = f"{blastp_result_dir}/blastp_result.xml"
                
                worker_logger.debug("Running DIAMOND to find homologues in FunBGCs")
                # Use parallel DIAMOND for better performance
                RunDIAMOND_parallel([combined_fasta], dmnd_db_SM_copy, blastp_result_dir, 1, 1, workers=2)
                
                worker_logger.debug("Adding homologue information")
                AddHomologue(blastp_result, df)
                
                # Duplication check
                dmnd_db_all = f"{dmnd_db_dir}/allCDS_{process_id}"
                
                worker_logger.debug("Creating DIAMOND database containing all proteins from the input")
                makeDIAMONDdb(fasta_file, dmnd_db_all)
                blastp_result2 = f"{blastp_result_dir}/blastp_result2.xml"
                
                worker_logger.debug("Running DIAMOND to find duplicated proteins")
                # Use parallel DIAMOND for better performance
                RunDIAMOND_parallel([combined_fasta], dmnd_db_all, blastp_result_dir, 2, 1, workers=2)
                
                worker_logger.debug("Checking duplications")
                DuplicationChecker(blastp_result2, min_identity, min_prot_len, df)
                
                # Remove small non-biosynthetic proteins from the dataframe
                worker_logger.debug("Removing small non-biosynthetic proteins")
                
                ToBeDeletedRows = []
                deletedProt = []
                for i in range(len(df)):
                    if df.at[i,"BP"] == 0 and df.at[i,"length"] < min_prot_len:
                        ToBeDeletedRows.append(i)
                        deletedProt.append(df.at[i,"locus_tag"])
                
                df_ = df.drop(index=ToBeDeletedRows).reset_index()
                original_df = df
                
                # Find clustered proteins
                worker_logger.debug("Finding clustered proteins")
                database = f"{current_dir}/data/DIAMOND/fungal_SM_proteins"
                ClusteredProteinFinder(extract_CDS_target_dir, database, IDdict, deletedProt, max_bgc_gap, df_, temp_dir_)
                
                # Extract BGCs
                BGC_dir = f"{temp_dir}/../BGCs"
                os.makedirs(BGC_dir, exist_ok=True)
                worker_logger.debug("Extracting BGCs")
                DefineBoundary(mode, GBK_dir, BGC_dir, gap_allowed, min_prot_len, fungus_name, df_, original_df, cluster_csv, all_cluster_csv, IDdict, GeneNumDict, MetabDict, temp_dir_, log_file)
                
                # Copy BGC gbk files
                if len(glob.glob(f"{BGC_dir}/*")) != 0:
                    fungus_name_ = fungus_name.replace(" ", "_")
                    worker_logger.debug("Copying BGC files")
                    all_BGC_dir = f"{results_dir}/all_clusters"
                    # Use lock to ensure safe directory operations
                    with directory_lock:
                        shutil.copytree(BGC_dir, f"{all_BGC_dir}/{fungus_name_}_BGCs", dirs_exist_ok=True)
                
                # Count BGCs
                BGCdf = pd.read_csv(cluster_csv)
                BGC_count = len(BGCdf)
                
            else:
                shutil.rmtree(output_dir)
                worker_logger.warning("No core/target protein was extracted")
                return (file_name, 0, output_dir, True, "No core/target protein was extracted")
        
        else:
            shutil.rmtree(output_dir)
            worker_logger.warning("The given input file does not contain any protein sequence.")
            return (file_name, 0, output_dir, True, "The given input file does not contain any protein sequence.")
        
        csv_dir = f"{output_dir}/CSVs"
        os.makedirs(csv_dir, exist_ok=True)
        analysis_csv = f"{csv_dir}/{file_name}_analysis.csv"
        analysis_csv2 = f"{csv_dir}/{file_name}_analysis_modified.csv"
        df.to_csv(analysis_csv)
        df_.to_csv(analysis_csv2)
        
        BGCdf = pd.read_csv(cluster_csv)
        BGCnum = len(BGCdf)
        if BGCnum > 0:
            worker_logger.debug("Creating HTML files")
            HTMLgenerator(cluster_csv, fungus_name_original, output_dir)
        else:
            shutil.rmtree(output_dir)
            worker_logger.warning("No BGC was extracted.")
            return (file_name, 0, output_dir, True, "No BGC was extracted.")
        
        # Delete temporary directory
        shutil.rmtree(temp_dir)
        
        worker_logger.info(f"{file_name}: BGC extraction finished. {BGCnum} BGCs were extracted.")
        
        return (file_name, BGCnum, output_dir, True, None)
    
    except Exception as e:
        # Log the error
        error_logger = logging.getLogger(f"error_{file_name}")
        error_logger.error(f"Error processing {file_name}: {str(e)}")
        return (file_name, 0, None, False, str(e))


def BGCeXtractor(gbk_dir,results_dir,mode,query,gap_allowed,max_bgc_gap,min_prot_len,num_of_genes_checked,min_identity,workers=None,max_memory=None,no_parallel=False):
    current_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Initialize logger first
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s:%(name)s - %(message)s')
    
    # Initialize memory monitor
    memory_monitor = MemoryMonitor(max_memory)
    
    # Log initial memory information
    mem_info = memory_monitor.get_memory_info()
    print(f"Memory info - Total: {mem_info['total']:.2f} GB, Available: {mem_info['available']:.2f} GB")
    if max_memory is not None:
        print(f"Memory limit set to: {max_memory} GB")
    
    # Determine parallel processing settings
    if no_parallel:
        use_parallel = False
        num_workers = 1
    else:
        use_parallel = True
        # Determine optimal number of workers based on system resources
        num_workers = get_optimal_worker_count(workers, max_memory)
    
    # Initialize resource locks for parallel processing
    manager = Manager()
    csv_lock = manager.Lock()  # For CSV file operations
    directory_lock = manager.Lock()  # For directory operations
    
    # Check if gbk_dir is a temporary directory created for FASTA conversion
    if os.path.basename(gbk_dir).startswith("temp_"):
        # Clean up temporary directory after processing
        cleanup_temp_dir = True
    else:
        cleanup_temp_dir = False

    """
    Create dictionaries
    """
    IDdict = {}
    GeneNumDict = {}
    MetabDict = {}
    protJSON = f"{current_dir}/data/proteins_indexed.json"
    MakeProtBGCidDict(IDdict,GeneNumDict,MetabDict,protJSON)

    """
    Create output directory
    """
    os.makedirs(results_dir,exist_ok=True)

    """
    Generate log file
    """
    log = logging.FileHandler(f'{results_dir}/{Path(gbk_dir).stem}.log')
    log.setFormatter(formatter)
    logger.addHandler(log)

    logger.info("BGC extraction started")

    """
    Collect GenBank files
    """
    input_gbk_files = glob.glob(f"{gbk_dir}/*")
    input_gbk_files.sort()

    # Filter out directories and only process files
    input_gbk_files = [f for f in input_gbk_files if os.path.isfile(f)]

    if len(input_gbk_files) == 0:
        logger.error("No GenBank files found in the input directory.")
        sys.exit("No GenBank files found in the input directory.")

    """
    Make a file to which the information of all clusters will be saved
    """
    columns = ["Fungus","BGC no.","Scaffold","Start position","End position","Locus tag start","Locus tag end",
                "Number of genes","Core enzymes","Similar BGC","Similarity score","Metabolite from similar BGC","Pfam domains"]

    logger.info(f"mode: {mode}")
    logger.info(f"query: {query}")
    logger.info(f"gap_allowed: {gap_allowed}")
    logger.info(f"max_bgc_gap: {max_bgc_gap}")
    logger.info(f"min_prot_len: {min_prot_len}")
    logger.info(f"num_of_genes_checked: {num_of_genes_checked}")
    logger.info(f"min_identity: {min_identity}")
    logger.info(f"parallel processing: {'enabled' if use_parallel else 'disabled'}")
    if use_parallel:
        logger.info(f"number of workers: {num_workers}")

    allBGCdf = pd.DataFrame(columns=columns)
    all_cluster_csv = f"{results_dir}/allBGCs.csv"
    # Use lock to ensure safe CSV file creation
    if use_parallel:
        with csv_lock:
            allBGCdf.to_csv(all_cluster_csv)
    else:
        allBGCdf.to_csv(all_cluster_csv)

    results_dir2 = f"{results_dir}/results"
    # Use lock to ensure safe directory creation
    if use_parallel:
        with directory_lock:
            os.makedirs(results_dir2, exist_ok=True)
    else:
        os.makedirs(results_dir2, exist_ok=True)

    all_BGC_dir = f"{results_dir}/all_clusters"
    # Use lock to ensure safe directory creation
    if use_parallel:
        with directory_lock:
            os.makedirs(all_BGC_dir, exist_ok=True)
    else:
        os.makedirs(all_BGC_dir, exist_ok=True)

    """
    Start analysis
    """
    if use_parallel and len(input_gbk_files) > 1:
        # Parallel processing with workload balancing
        logger.info(f"Starting parallel processing of {len(input_gbk_files)} files using {num_workers} workers")
        
        # Estimate file sizes for workload balancing
        file_sizes = []
        for file in input_gbk_files:
            try:
                file_size = os.path.getsize(file)
                file_sizes.append((file, file_size))
            except:
                # Default to 1MB if size cannot be determined
                file_sizes.append((file, 1024*1024))
        
        # Sort files by size (largest first) for better workload distribution
        file_sizes.sort(key=lambda x: x[1], reverse=True)
        
        # Prepare arguments for worker processes with workload balancing and resource locks
        worker_args = []
        for file, _ in file_sizes:
            worker_args.append((file, gbk_dir, results_dir2, mode, query, gap_allowed,
                               max_bgc_gap, min_prot_len, num_of_genes_checked, min_identity,
                               current_dir, IDdict, GeneNumDict, MetabDict, all_cluster_csv, csv_lock, directory_lock))
        
        # Process files in parallel with progress tracking and error handling
        # Using chunksize for better workload distribution
        chunksize = max(1, len(input_gbk_files) // (num_workers * 2))
        with Pool(processes=num_workers) as pool:
            try:
                results = list(tqdm(pool.imap(process_single_file_worker, worker_args, chunksize=chunksize),
                                    total=len(input_gbk_files),
                                    desc="Processing GenBank files"))
            except Exception as e:
                logger.error(f"Error in parallel processing: {str(e)}")
                # Fall back to sequential processing if parallel processing fails
                logger.warning("Falling back to sequential processing due to error")
                use_parallel = False
                results = []
                for args in worker_args:
                    result = retry_operation(process_single_file_worker, (args,), max_retries=2)
                    if result is None:
                        file_name = args[0]
                        results.append((file_name, 0, None, False, "Processing failed after retries"))
                    else:
                        results.append(result)
        
        # Process results
        total_BGCs = 0
        successful_files = 0
        failed_files = 0
        
        for file_name, BGC_count, output_dir, success, error_message in results:
            if success:
                if BGC_count > 0:
                    successful_files += 1
                    total_BGCs += BGC_count
                    logger.info(f"{file_name}: Successfully processed. {BGC_count} BGCs extracted.")
                else:
                    logger.warning(f"{file_name}: Processed but no BGCs extracted.")
            else:
                failed_files += 1
                logger.error(f"{file_name}: Processing failed. Error: {error_message}")
        
        logger.info(f"Parallel processing completed. {successful_files}/{len(input_gbk_files)} files processed successfully.")
        logger.info(f"Total BGCs extracted: {total_BGCs}")
        if failed_files > 0:
            logger.warning(f"{failed_files} files failed to process.")
    
    else:
        # Sequential processing (original approach)
        logger.info("Starting sequential processing of GenBank files")
        
        for file in input_gbk_files:
            time_start = time.time()

            file_name = Path(file).stem
            output_dir = f"{results_dir2}/{file_name}_results"
            temp_dir = f"{output_dir}/temp"
            temp_dir_ = f"{temp_dir}/temp"
            os.makedirs(output_dir, exist_ok=True)
            os.makedirs(temp_dir, exist_ok=True)
            df = pd.DataFrame()

            print(f"{file_name}: BGC extraction started")
            logger.info(f"{file_name}: BGC extraction started")


            """
            Make a file to which the information of all clusters will be saved
            """
            columns = ["BGC no.","Scaffold","Start position","End position","Locus tag start","Locus tag end",
                        "Number of genes","Core enzymes","Similar BGC","Similarity score","Metabolite from similar BGC","Pfam domains"]

            logger.info(f"{mode} mode was selected")

            BGCdf = pd.DataFrame(columns=columns)
            cluster_csv = f"{output_dir}/BGCs.csv"
            BGCdf.to_csv(cluster_csv)

            """
            Check the input file. If locus tags are not given for CDSs, they are added as PROTEINXXXXXX.
            """
            gbk_file = SeqHandling(file,df)
            gbk_file.AddLocusTag(temp_dir)
            modified_gbk_file = glob.glob(f"{temp_dir}/*.gbk")[0]
            modified_gbk = SeqHandling(modified_gbk_file,df)

            logger.debug("Processed the input GenBank file")

            """
            Make necessary directories
            """
            CDS_dir = f"{temp_dir}/extracted_CDS"
            GBK_dir = f"{temp_dir}/scaffolds"

            os.makedirs(CDS_dir,exist_ok=True)
            os.makedirs(GBK_dir,exist_ok=True)

            
            """
            Obtain the fungus's name
            """ 
            fungus_name = modified_gbk.name() # Special symbols that might cause errors removed
            fungus_name_original = modified_gbk.name_o() # Original fungus's name

            logger.debug(f"Original fungus's name: {fungus_name_original}")
            logger.debug(f"Modified version of fungus's name: {fungus_name}")


            """
            Extract all CDSs as a single fasta file and each scaffold as a gbk file
            """
            modified_gbk.CDS_extractor(CDS_dir)
            modified_gbk.GBK_extractor(GBK_dir)
            fasta_file = f"{CDS_dir}/{file_name}.fasta"

            logger.debug("All protein sequences were extracted")
            logger.debug("All scaffolds were extracted")

            if os.stat(fasta_file).st_size != 0:
                        
                dmnd_db_dir = f"{temp_dir}/BLASTdb"
                os.makedirs(dmnd_db_dir,exist_ok=True)

                if mode == "all":
                    """
                    hmmscan to detect core enzymes
                    """
                    hmmscan_output_dir = f"{temp_dir}/hmmscan"
                    os.makedirs(hmmscan_output_dir,exist_ok=True)
                    hmmscan_result = f"{hmmscan_output_dir}/hmm_core.txt"
                    hmmscan_core_database = f"{current_dir}/data/hmm/core/core_UstY.hmm"
                    
                    logger.debug("Running hmmscan to detect core enzymes")
                    runHMMscan(fasta_file,hmmscan_result,hmmscan_core_database,"1e-5")
                    
                    logger.debug("Finding core enzymes")
                    Core_extractor(hmmscan_result,df,temp_dir_)

                if mode == "target":
                    """
                    Create DIAMOND database
                    """
                    target_dmnd_db = f"{dmnd_db_dir}/target"
                    logger.debug("Creating DIAMOND database using the input fasta file")
                    makeDIAMONDdb(query,target_dmnd_db)

                    """
                    Perform blastp to find target homologues
                    """   
                    blastp_result_dir = f"{temp_dir}/blastp_results"
                    os.makedirs(blastp_result_dir,exist_ok=True)

                    target_blastp_result = f"{blastp_result_dir}/target_blastp_result.xml"
                    logger.debug("Running DIAMOND to find target proteins")
                    RunDIAMOND(fasta_file,target_dmnd_db,target_blastp_result,1,1)

                    """
                    Add target homologue information to the csv file
                    """
                    logger.debug("Adding the information of the closest homologue of each target")
                    AddTargetHomologue(target_blastp_result,df)

                if mode == "pfam":
                    """
                    Extract the target Pfam HMM file
                    """
                    funSMhmm = f"{current_dir}/data/hmm/fungalSM/fungalSM.hmm"
                    hmm_dir = f"{temp_dir}/hmm"
                    os.makedirs(hmm_dir,exist_ok=True)
                    target_hmm = f"{hmm_dir}/{query}.hmm"

                    logger.debug("Fetching the required HMM file")
                    runHMMfetch(funSMhmm,query,target_hmm)

                    if os.stat(target_hmm).st_size == 0:
                        logger.error("The given protein family name was not found in the database.")
                        sys.exit("The given protein family name was not found in the database.")
                    
                    logger.debug("Running hmmpress")
                    runHMMpress(target_hmm)

                    """
                    hmmscan to detect target proteins
                    """
                    hmmscan_output_dir = f"{temp_dir}/hmmscan"
                    os.makedirs(hmmscan_output_dir,exist_ok=True)
                    hmmscan_result = f"{hmmscan_output_dir}/hmm_target.txt"
                    hmmscan_target_database = f"{temp_dir}/hmm/{query}.hmm"

                    logger.debug("Running hmmscan to detect target proteins")
                    runHMMscan(fasta_file,hmmscan_result,hmmscan_target_database,"1e-5")

                    logger.debug("Adding target information")
                    AddTargetPfam(hmmscan_result,df)

                """
                Extract protein sequences around core/target proteins
                """
                extract_CDS_target_dir = f"{temp_dir}/extracted_CDS_target"
                os.makedirs(extract_CDS_target_dir,exist_ok=True)

                logger.debug("Extracting protein sequences around core/target proteins")
                ExtractCDS4Check(mode,num_of_genes_checked,extract_CDS_target_dir,df)

                """
                Combine all extracted fasta files as a single file
                """
                fasta_files = glob.glob(f"{extract_CDS_target_dir}/*.fasta")
                if len(fasta_files) != 0:
                    logger.debug("Combining all fasta files")
                    combineFASTA(extract_CDS_target_dir,temp_dir)
                    combined_fasta = f"{temp_dir}/combined.fasta"

                    if mode == "target" or mode == "pfam":
                        """
                        hmmscan to detect core enzymes
                        """
                        hmmscan_output_dir = f"{temp_dir}/hmmscan"
                        os.makedirs(hmmscan_output_dir,exist_ok=True)
                        hmmscan_result = f"{hmmscan_output_dir}/hmm_core.txt"
                        hmmscan_core_database = f"{current_dir}/data/hmm/core/core.hmm"

                        logger.debug("Running hmmscan to detect core enzymes")
                        runHMMscan(combined_fasta,hmmscan_result,hmmscan_core_database,"1e-5")

                        logger.debug("Finding core enzymes")
                        Core_extractor(hmmscan_result,df,temp_dir_)

                    logger.debug("Finding potential RiPP precursor proteins")
                    aa_len, unique_aa, min_repeat, min_match = 10, 3, 3, 9
                    RiPPppFinder(aa_len,unique_aa,min_repeat,min_match,combined_fasta,df)

                    """
                    Perform hmmscan using the fungal SM protein domain database
                    """        
                    hmmscan_result_ = f"{hmmscan_output_dir}/hmmscan_result.txt"
                    hmmscan_SM_database = f"{current_dir}/data/hmm/fungalSM/fungalSM.hmm"

                    logger.debug("Running hmmscan to detect potential fungal SM proteins")
                    runHMMscan(combined_fasta,hmmscan_result_,hmmscan_SM_database,"1e-5")

                    """
                    Add Pfam domain information to the csv file
                    """
                    SMhmm = f"{current_dir}/data/hmm/fungalSM/fungalSM.hmm"

                    logger.debug("Adding Pfam information")
                    AddPfam(hmmscan_result_,df,SMhmm)

                    """
                    Find the closest homologue in the fungal biosynthetic proteins
                    """
                    dmnd_db_SM = f"{current_dir}/data/DIAMOND/fungal_SM_proteins"
                    blastp_result_dir = f"{temp_dir}/blastp_results"
                    os.makedirs(blastp_result_dir,exist_ok=True)
                    blastp_result = f"{blastp_result_dir}/blastp_result.xml"

                    logger.debug("Running DIAMOND to find homologues in FunBGCs")
                    RunDIAMOND(combined_fasta,dmnd_db_SM,blastp_result,1,1)

                    logger.debug("Adding homologue information")
                    AddHomologue(blastp_result,df)

                    """
                    Duplication check
                    """
                    dmnd_db_all = f"{dmnd_db_dir}/allCDS"

                    logger.debug("Creating DIAMOND database containing all proteins from the input")
                    makeDIAMONDdb(fasta_file,dmnd_db_all)      
                    blastp_result2 = f"{blastp_result_dir}/blastp_result2.xml"

                    logger.debug("Running DIAMOND to find duplicated proteins")
                    RunDIAMOND(combined_fasta,dmnd_db_all,blastp_result2,2,1)

                    logger.debug("Checking duplications")
                    DuplicationChecker(blastp_result2,min_identity,min_prot_len,df)

                    """
                    Remove small non-biosynthetic proteins from the dataframe
                    """        
                    logger.debug("Removing small non-biosynthetic proteins")
                    
                    ToBeDeletedRows = []
                    deletedProt = []
                    for i in range(len(df)):
                        if df.at[i,"BP"] == 0 and df.at[i,"length"] < min_prot_len:
                            ToBeDeletedRows.append(i)
                            deletedProt.append(df.at[i,"locus_tag"])

                    df_ = df.drop(index=ToBeDeletedRows).reset_index()
                    original_df = df

                    """
                    Find clustered proteins
                    """
                    logger.debug("Finding clustered proteins")
                    database = f"{current_dir}/data/DIAMOND/fungal_SM_proteins"
                    ClusteredProteinFinder(extract_CDS_target_dir,database,IDdict,deletedProt,max_bgc_gap,df_,temp_dir_)

                    """
                    Extract BGCs
                    """
                    BGC_dir = f"{temp_dir}/../BGCs"
                    os.makedirs(BGC_dir,exist_ok=True)
                    logger.debug("Extracting BGCs")
                    DefineBoundary(mode,GBK_dir,BGC_dir,gap_allowed,min_prot_len,fungus_name,df_,original_df,cluster_csv,all_cluster_csv,IDdict,GeneNumDict,MetabDict,temp_dir_,log)

                    """
                    Copy BGC gbk files
                    """
                    if len(glob.glob(f"{BGC_dir}/*")) != 0:
                        fungus_name_ = fungus_name.replace(" ","_")
                        logger.debug("Copying BGC files")
                        shutil.copytree(BGC_dir,f"{all_BGC_dir}/{fungus_name_}_BGCs",dirs_exist_ok=True)

                else:
                    shutil.rmtree(output_dir)
                    logger.warning("No core/target protein was extracted")
                    print("No core/target protein was extracted")
                    continue

            else:
                shutil.rmtree(output_dir)
                logger.warning("The given input file does not contain any protein sequence.")
                print("The given input file does not contain any protein sequence.")
                continue

            csv_dir = f"{output_dir}/CSVs"
            os.makedirs(csv_dir,exist_ok=True)
            analysis_csv = f"{csv_dir}/{file_name}_analysis.csv"
            analysis_csv2 = f"{csv_dir}/{file_name}_analysis_modified.csv"
            df.to_csv(analysis_csv)
            df_.to_csv(analysis_csv2)

            BGCdf = pd.read_csv(cluster_csv)
            BGCnum = len(BGCdf)
            if BGCnum > 0:
                logger.debug("Creating HTML files")
                HTMLgenerator(cluster_csv,fungus_name_original,output_dir)
            else:
                shutil.rmtree(output_dir)
                logger.warning("No BGC was extracted.")
                print("No BGC was extracted.")
                continue

            time_end = time.time()
            total_time = "{:.2f}".format((time_end - time_start)/60)

            """
            Delete directories
            """
            shutil.rmtree(temp_dir)


            """
            Message for BGC extraction completion
            """    
            if BGCnum == 1:
                message = f"{BGCnum} BGC was"
            else:
                message = f"{BGCnum} BGCs were"

            logger.info(f"{file_name}: BGC extraction finished in {total_time} min. {message} extracted.")
            print(f"{file_name}: BGC extraction finished in {total_time} min. {message} extracted.")


    logger.info("All BGC extractions completed")
    
    # Copy GenBank file to output directory if it was created for FASTA conversion
    if cleanup_temp_dir:
        # Find the GenBank file in the temporary directory
        gbk_files = glob.glob(f"{gbk_dir}/*.gbk")
        if gbk_files:
            # Copy the GenBank file to the output directory
            for gbk_file in gbk_files:
                shutil.copy2(gbk_file, results_dir)
                logger.info(f"Copied GenBank file to output directory: {os.path.join(results_dir, os.path.basename(gbk_file))}")
        
        # Clean up temporary directory
        shutil.rmtree(gbk_dir)
        logger.debug(f"Cleaned up temporary directory: {gbk_dir}")