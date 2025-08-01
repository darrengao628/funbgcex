#!/usr/bin/env python3
"""
CSV to JSON Converter for FunBGCeX Protein Database

This script converts the protein database from CSV format to JSON format
for improved performance and easier data manipulation in Python.
"""

import pandas as pd
import json
import os
import time
from pathlib import Path

def convert_csv_to_json(csv_file_path, json_file_path, indexed_json_file_path=None):
    """
    Convert CSV file to JSON format.
    
    Args:
        csv_file_path (str): Path to the input CSV file
        json_file_path (str): Path to the output JSON file
        indexed_json_file_path (str, optional): Path to the output indexed JSON file
    """
    start_time = time.time()
    
    print(f"Reading CSV file: {csv_file_path}")
    # Read the CSV file using pandas for efficient processing
    df = pd.read_csv(csv_file_path)
    
    # Convert DataFrame to list of dictionaries (JSON-like structure)
    print("Converting to JSON format...")
    data = df.to_dict(orient='records')
    
    # Write the standard JSON file
    print(f"Writing JSON file: {json_file_path}")
    with open(json_file_path, 'w') as json_file:
        json.dump(data, json_file, indent=2)
    
    # Create and write the indexed JSON file if requested
    if indexed_json_file_path:
        print(f"Creating indexed JSON file: {indexed_json_file_path}")
        
        # Create a dictionary indexed by protein ID for faster lookups
        indexed_data = {}
        for record in data:
            protein_id = record.get('protein ID', '')
            if protein_id:
                indexed_data[protein_id] = record
        
        # Write the indexed JSON file
        with open(indexed_json_file_path, 'w') as indexed_json_file:
            json.dump(indexed_data, indexed_json_file, indent=2)
    
    end_time = time.time()
    print(f"Conversion completed in {end_time - start_time:.2f} seconds")
    print(f"Total records processed: {len(data)}")

def main():
    # Determine paths relative to this script
    script_dir = Path(__file__).parent
    data_dir = script_dir.parent / 'data'
    
    # Input and output file paths
    csv_file_path = data_dir / 'proteins.csv'
    json_file_path = data_dir / 'proteins.json'
    indexed_json_file_path = data_dir / 'proteins_indexed.json'
    
    # Check if input file exists
    if not csv_file_path.exists():
        print(f"Error: Input file not found: {csv_file_path}")
        return
    
    # Convert CSV to JSON
    convert_csv_to_json(
        csv_file_path=csv_file_path,
        json_file_path=json_file_path,
        indexed_json_file_path=indexed_json_file_path
    )
    
    print("\nConversion completed successfully!")
    print(f"Standard JSON file: {json_file_path}")
    print(f"Indexed JSON file: {indexed_json_file_path}")

if __name__ == "__main__":
    main()