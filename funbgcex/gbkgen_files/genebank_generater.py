#!/usr/bin/env python3
 
import argparse
import os
import subprocess
import sys
from pathlib import Path
import logging
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed
 
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
 
# Import necessary modules from antiSMASH codebase
from . import gff_parser
 
# Custom exception class to replace AntismashInputError
class AntismashInputError(Exception):
    """Exception raised for invalid input in GFF parsing."""
    pass
 
def setup_logging():
    """Set up logging configuration."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    return logging.getLogger('fasta_to_genbank')
 
def run_augustus_single_sequence(args):
    """Run Augustus on a single sequence (helper function for multiprocessing).
    
    Args:
        args: Tuple of (sequence_record, output_dir, species, index)
        
    Returns:
        Path to the generated GFF3 file for this sequence
    """
    record, output_dir, species, index = args
    logger = logging.getLogger('augustus')
    
    # Create individual FASTA file for this sequence
    fasta_file = os.path.join(output_dir, f"sequence_{index}.fasta")
    with open(fasta_file, 'w') as f:
        SeqIO.write(record, f, 'fasta')
    
    # Set output GFF file path
    gff_file = os.path.join(output_dir, f"sequence_{index}.gff3")
    
    # Augustus command
    cmd = [
        "augustus",
        f"--species={species}",
        "--gff3=on",
        "--UTR=off",
        "--protein=on",  # Add protein output
        "--outfile=" + gff_file,
        fasta_file
    ]
    
    logger.info(f"Running Augustus on sequence {index}: {record.id}")
    
    try:
        process = subprocess.run(
            cmd,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        logger.info(f"Augustus finished successfully for sequence {index}")
        return gff_file
    except subprocess.CalledProcessError as e:
        logger.error(f"Augustus failed for sequence {index}: {e.stderr}")
        raise RuntimeError(f"Gene prediction with Augustus failed for sequence {index}") from e
 
def run_augustus(fasta_file, output_dir, species="Aspergillus fumigatus", cpu_cores=None):
    """Run Augustus to predict genes and generate GFF3 output with multiprocessing.
    
    Args:
        fasta_file: Path to input FASTA file
        output_dir: Directory for Augustus output
        species: Augustus species model to use (default: generic)
        cpu_cores: Number of CPU cores to use (default: all available)
        
    Returns:
        Path to the combined GFF3 file
    """
    logger = logging.getLogger('augustus')
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Load FASTA records
    records = list(SeqIO.parse(fasta_file, "fasta"))
    logger.info(f"Processing {len(records)} sequences with Augustus")
    
    if len(records) == 1 or cpu_cores == 1:
        # Single sequence or single core - use original method
        gff_file = os.path.join(output_dir, f"{Path(fasta_file).stem}.gff3")
        cmd = [
            "augustus",
            f"--species={species}",
            "--gff3=on",
            "--UTR=off",
            "--protein=on",  # Add protein output
            "--outfile=" + gff_file,
            fasta_file
        ]
        
        logger.info(f"Running Augustus: {' '.join(cmd)}")
        try:
            process = subprocess.run(
                cmd,
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )
            logger.info("Augustus finished successfully")
            return gff_file
        except subprocess.CalledProcessError as e:
            logger.error(f"Augustus failed: {e.stderr}")
            raise RuntimeError("Gene prediction with Augustus failed") from e
    
    # Multiprocessing for multiple sequences
    if cpu_cores is None:
        cpu_cores = min(multiprocessing.cpu_count(), len(records))
    
    logger.info(f"Using {cpu_cores} CPU cores for multiprocessing")
    
    # Prepare arguments for multiprocessing
    args_list = [(record, output_dir, species, i) for i, record in enumerate(records)]
    
    # Run Augustus in parallel
    gff_files = []
    with ProcessPoolExecutor(max_workers=cpu_cores) as executor:
        future_to_gff = {executor.submit(run_augustus_single_sequence, args): args[3]
                        for args in args_list}
        
        for future in as_completed(future_to_gff):
            index = future_to_gff[future]
            try:
                gff_file = future.result()
                gff_files.append(gff_file)
            except Exception as e:
                logger.error(f"Error processing sequence {index}: {e}")
                raise
    
    # Combine all GFF files into one
    combined_gff = os.path.join(output_dir, f"{Path(fasta_file).stem}.gff3")
    with open(combined_gff, 'w') as outfile:
        for gff_file in sorted(gff_files):
            with open(gff_file, 'r') as infile:
                outfile.write(infile.read())
                outfile.write("\n")
    
    # Clean up temporary files
    for gff_file in gff_files:
        try:
            os.remove(gff_file)
            os.remove(gff_file.replace('.gff3', '.fasta'))
        except:
            pass
    
    logger.info(f"Combined GFF file created: {combined_gff}")
    return combined_gff
 
def load_fasta_records(fasta_file):
    """Load FASTA file into BioPython records.
    
    Args:
        fasta_file: Path to FASTA file
        
    Returns:
        List of BioPython SeqRecord objects
    """
    logger = logging.getLogger('fasta_loader')
    logger.info(f"Loading FASTA file: {fasta_file}")
    
    try:
        records = list(SeqIO.parse(fasta_file, "fasta"))
        # Set molecule_type annotation required for GenBank
        for record in records:
            record.annotations["molecule_type"] = "DNA"
        logger.info(f"Loaded {len(records)} records from FASTA file")
        return records
    except Exception as e:
        logger.error(f"Failed to load FASTA: {e}")
        raise
 
def convert_to_genbank(fasta_records, gff_file, output_file):
    """Convert FASTA records and GFF annotations to GenBank format.
    
    Args:
        fasta_records: List of BioPython SeqRecord objects
        gff_file: Path to GFF3 file from Augustus
        output_file: Path to write the output GenBank file
    """
    logger = logging.getLogger('genbank_converter')
    
    try:
        # Check if GFF is suitable for the FASTA records
        logger.info("Checking GFF suitability")
        gff_parser.check_gff_suitability(gff_file, fasta_records)
        
        # Add features from GFF to records
        logger.info("Adding features from GFF to records")
        gff_parser.update_records(gff_file, fasta_records)
        
        
        # Write records to GenBank format
        logger.info(f"Writing GenBank output to {output_file}")
        with open(output_file, 'w') as handle:
            SeqIO.write(fasta_records, handle, 'genbank')
            
        logger.info("GenBank conversion completed successfully")
    
    except Exception as e:
        logger.error(f"GenBank conversion failed: {str(e)}")
        raise
 
def main(args=None):
    """Main function to handle the conversion process.
    
    Args:
        args: List of command line arguments (without script name).
             If None, uses sys.argv[1:]
    """
    parser = argparse.ArgumentParser(description="Convert DNA FASTA to GenBank with Augustus gene prediction")
    
    parser.add_argument("input", help="Input DNA FASTA file")
    parser.add_argument("-o", "--output", help="Output GenBank file (optional, defaults to input filename with .gbk extension)")
    parser.add_argument("-s", "--species", default="aspergillus_fumigatus",
                       help="Augustus species model (default: generic)")
    parser.add_argument("-w", "--workdir", default="./augustus_output",
                       help="Working directory for Augustus (default: ./augustus_output)")
    parser.add_argument("-g", "--gff", help="Optional pre-existing GFF3 file to use instead of running Augustus")
    parser.add_argument("-c", "--cpu", type=int, default=None,
                       help="Number of CPU cores to use for Augustus (default: all available)")
    
    args = parser.parse_args(args)
    
    # Set default output file name if not specified
    if not args.output:
        args.output = os.path.splitext(args.input)[0] + ".gbk"
    
    # Check for corresponding GFF file if not explicitly provided
    if not args.gff:
        input_base = os.path.splitext(args.input)[0]
        possible_gff = input_base + ".gff"
        possible_gff3 = input_base + ".gff3"
        
        if os.path.exists(possible_gff):
            args.gff = possible_gff
            logger = setup_logging()
            logger.info(f"Found corresponding GFF file: {possible_gff}")
        elif os.path.exists(possible_gff3):
            args.gff = possible_gff3
            logger = setup_logging()
            logger.info(f"Found corresponding GFF3 file: {possible_gff3}")
    
    logger = setup_logging()
    logger.info(f"Starting conversion of {args.input} to {args.output}")
    
    try:
        # Load FASTA file as BioPython records
        fasta_records = load_fasta_records(args.input)
        
        # Use provided GFF file or run Augustus to generate one
        if args.gff:
            if not os.path.exists(args.gff):
                raise FileNotFoundError(f"GFF file not found: {args.gff}")
            gff_file = args.gff
            logger.info(f"Using provided GFF file: {gff_file}")
        else:
            gff_file = run_augustus(args.input, args.workdir, args.species, args.cpu)
            logger.info(f"Generated GFF file from Augustus: {gff_file}")
        
        # Convert to GenBank
        convert_to_genbank(fasta_records, gff_file, args.output)
        
        logger.info(f"Successfully created GenBank file: {args.output}")
        return 0
        
    except Exception as e:
        logger.error(f"Conversion failed: {e}")
        return 1
 
if __name__ == "__main__":
    sys.exit(main())