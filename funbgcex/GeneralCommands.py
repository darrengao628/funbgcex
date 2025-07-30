#!/usr/bin/env python3

import subprocess
import os
import tempfile
import glob
from multiprocessing import Pool, cpu_count
from functools import partial


def runHMMscan(input_file,output_file,database,Evalue):
    subprocess.run(["hmmscan","--domtblout",output_file,"-E",Evalue,"--domE",
                    Evalue,database,input_file],stdout=subprocess.DEVNULL)

def runHMMfetch(HMMdatabase,target_name,output_file):
    subprocess.run(["hmmfetch","-o",output_file,HMMdatabase,target_name],stdout=subprocess.DEVNULL)

def runHMMpress(input_file):
    subprocess.run(["hmmpress",input_file],stdout=subprocess.DEVNULL)

def makeDIAMONDdb(fasta,database):
    subprocess.run(["diamond","makedb","--in",fasta,"--db",database],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)

def RunDIAMOND(fasta,database,output_file,max_target_seqs,max_hsps):
    max_target_seqs = str(max_target_seqs)
    max_hsps = str(max_hsps)
    subprocess.run(["diamond","blastp","--query",fasta,"--db",database,
                    "--out",output_file,"--outfmt","5","--max-target-seqs",max_target_seqs,
                    "--max-hsps",max_hsps],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)

def RunDIAMOND_parallel(fasta_list, database, output_dir, max_target_seqs, max_hsps, workers=None):
    """
    Run DIAMOND searches in parallel for multiple FASTA files.
    
    Args:
        fasta_list: List of paths to FASTA files to search
        database: Path to DIAMOND database
        output_dir: Directory to store output files
        max_target_seqs: Maximum number of target sequences to report
        max_hsps: Maximum number of high-scoring segment pairs to report
        workers: Number of worker processes (default: auto-detect)
    
    Returns:
        List of paths to output files
    """
    if workers is None:
        workers = max(1, cpu_count() - 1)
    
    max_target_seqs = str(max_target_seqs)
    max_hsps = str(max_hsps)
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Prepare arguments for each worker
    worker_args = []
    output_files = []
    
    for fasta in fasta_list:
        fasta_name = os.path.splitext(os.path.basename(fasta))[0]
        output_file = os.path.join(output_dir, f"{fasta_name}_blastp_result.xml")
        output_files.append(output_file)
        worker_args.append((fasta, database, output_file, max_target_seqs, max_hsps))
    
    # Define worker function
    def _run_diamond_worker(args):
        fasta, db, out, max_targets, max_h = args
        subprocess.run(["diamond", "blastp", "--query", fasta, "--db", db,
                        "--out", out, "--outfmt", "5", "--max-target-seqs", max_targets,
                        "--max-hsps", max_h], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return out
    
    # Process files in parallel
    with Pool(processes=workers) as pool:
        results = pool.map(_run_diamond_worker, worker_args)
    
    return output_files

def runHMMscan_parallel(input_files, database, output_dir, Evalue, workers=None):
    """
    Run hmmscan searches in parallel for multiple input files.
    
    Args:
        input_files: List of paths to input files to search
        database: Path to HMM database
        output_dir: Directory to store output files
        Evalue: E-value threshold for reporting hits
        workers: Number of worker processes (default: auto-detect)
    
    Returns:
        List of paths to output files
    """
    if workers is None:
        workers = max(1, cpu_count() - 1)
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Prepare arguments for each worker
    worker_args = []
    output_files = []
    
    for input_file in input_files:
        input_name = os.path.splitext(os.path.basename(input_file))[0]
        output_file = os.path.join(output_dir, f"{input_name}_hmmscan_result.txt")
        output_files.append(output_file)
        worker_args.append((input_file, output_file, database, Evalue))
    
    # Define worker function
    def _run_hmmscan_worker(args):
        input_f, out, db, e_val = args
        subprocess.run(["hmmscan", "--domtblout", out, "-E", e_val, "--domE",
                        e_val, db, input_f], stdout=subprocess.DEVNULL)
        return out
    
    # Process files in parallel
    with Pool(processes=workers) as pool:
        results = pool.map(_run_hmmscan_worker, worker_args)
    
    return output_files
