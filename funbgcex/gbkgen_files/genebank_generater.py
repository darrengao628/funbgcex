#!/usr/bin/env python3

"""
Generator for synthetic GenBank files for testing purposes.
"""

import random
import string
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO

def generate_random_dna(length):
    """
    Generate a random DNA sequence.
    
    Args:
        length: Length of the DNA sequence
        
    Returns:
        Random DNA sequence as a string
    """
    return ''.join(random.choice('ATCG') for _ in range(length))

def generate_random_protein(length):
    """
    Generate a random protein sequence.
    
    Args:
        length: Length of the protein sequence
        
    Returns:
        Random protein sequence as a string
    """
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    return ''.join(random.choice(amino_acids) for _ in range(length))

def generate_test_genbank(output_file, num_genes=100, scaffold_length=50000):
    """
    Generate a synthetic GenBank file for testing.
    
    Args:
        output_file: Path to output GenBank file
        num_genes: Number of genes to generate
        scaffold_length: Length of the scaffold
    """
    # Generate random scaffold sequence
    scaffold_seq = generate_random_dna(scaffold_length)
    
    # Create SeqRecord
    record = SeqRecord(
        Seq(scaffold_seq),
        id=f"test_scaffold_1",
        name=f"TestFungus",
        description="Synthetic test genome for FunBGCeX performance testing",
        annotations={"molecule_type": "DNA"}
    )
    
    # Generate genes
    gene_length = scaffold_length // (num_genes + 1)  # Leave space between genes
    gene_positions = []
    
    for i in range(num_genes):
        # Calculate gene position with some spacing
        start = i * gene_length + random.randint(100, 500)
        end = start + random.randint(300, 3000)  # Gene length between 300-3000 bp
        
        # Ensure we don't go beyond scaffold length
        if end >= scaffold_length:
            break
            
        gene_positions.append((start, end))
        
        # Generate CDS feature
        feature = SeqFeature(
            FeatureLocation(start, end, strand=1),
            type="CDS",
            qualifiers={
                "locus_tag": f"PROTEIN{i+1:06d}",
                "product": f"hypothetical protein {i+1}",
                "translation": generate_random_protein((end - start) // 3)
            }
        )
        
        record.features.append(feature)
    
    # Write to file
    with open(output_file, "w") as output_handle:
        SeqIO.write(record, output_handle, "genbank")

def main(args):
    """
    Main function for genebank generator.
    
    Args:
        args: List of command line arguments
    """
    import argparse
    
    parser = argparse.ArgumentParser(description="Generate synthetic GenBank files")
    parser.add_argument("input_file", help="Input FASTA file")
    parser.add_argument("-o", "--output", help="Output GenBank file", required=True)
    parser.add_argument("--num-genes", type=int, default=100, help="Number of genes to generate")
    parser.add_argument("--scaffold-length", type=int, default=50000, help="Length of the scaffold")
    
    parsed_args = parser.parse_args(args)
    
    # For now, just create a synthetic GenBank file
    # In a real implementation, this would convert the FASTA file
    generate_test_genbank(parsed_args.output, parsed_args.num_genes, parsed_args.scaffold_length)
    print(f"Generated test GenBank file: {parsed_args.output}")

if __name__ == "__main__":
    # Test the generator
    generate_test_genbank("test_genbank.gbk", num_genes=10)
    print("Generated test GenBank file: test_genbank.gbk")