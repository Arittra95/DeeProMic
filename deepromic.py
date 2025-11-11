#!/usr/bin/env python
# Welcome to DeeproMic! 
# File name         :deepromic.py
# description       :Classification of druggable targets 
# author            :Arittra Bhattacharjee; email: arittra.bioinfo@gmail.com 
# date of creation  :17 December 2023 to 19 December 2023 
# version           :1.0
# command           :python deepromic.py -i <path_to_input_fasta> -o <path_to_output_directory>
# required files    :human.fasta, essential.fasta and profeatx (executable file)  
# conda env create -f environment.yml -n deepromic
#**************************************************************************************************

import os
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"
import tensorflow as tf
tf.config.list_physical_devices("GPU")
import argparse
import subprocess
import pandas as pd
import numpy as np
import joblib
from tensorflow.keras.models import load_model
from Bio import SeqIO
import warnings
warnings.filterwarnings("ignore")

print("Welcome to DeeProMic!")

def main():
    parser = argparse.ArgumentParser(description='DeeProMic: Therapeutic Protein Classifier')
    parser.add_argument('-i', '--input', required=True, help='Input protein FASTA file')
    parser.add_argument('-t', '--threshold', type=float, default=0.5, help='Probability threshold for Class 1')
    parser.add_argument('-o', '--output', default='output', help='Output directory')
    
    args = parser.parse_args()
    input_fasta = args.input
    threshold = args.threshold
    output_directory = args.output

    # Create output directory
    os.makedirs(output_directory, exist_ok=True)
    print(f"Output directory: {output_directory}")

    # Clean FASTA headers
    base_name = os.path.splitext(os.path.basename(input_fasta))[0]
    modified_input = os.path.join(output_directory, f"{base_name}_modified.fasta")
    print("Cleaning FASTA headers...")
    modify_cmd = f"awk '/^>/{{sub(/ .*/, \"\"); print; next}} {{if(NF) print}}' {input_fasta} > {modified_input}"
    subprocess.run(modify_cmd, shell=True, check=True)

    # Generate DDE features
    print("Converting sequences to DDE features...")
    dde_output = os.path.join(output_directory, "dde.csv")
    profeat_cmd = f"./profeatx -i {modified_input} -o {dde_output} -e DDE"
    subprocess.run(profeat_cmd, shell=True, check=True)

    # Load model and scaler
    print("Loading model and scaler...")
    model = load_model("lstm_trained_model.h5")
    scaler = joblib.load("lstm_scaler_model.joblib")

    # Read and process DDE features
    df = pd.read_csv(dde_output, sep='\t')
    df2 = pd.read_csv(dde_output, sep='\t', header=0, index_col=0)
    df_features = df.iloc[:, 1:-1]

    # Make predictions
    X = df_features.values
    X_scaled = scaler.transform(X)
    X_reshaped = np.reshape(X_scaled, (X_scaled.shape[0], 1, X_scaled.shape[1]))
    predictions = model.predict(X_reshaped)

    # Create results dataframe
    rounded_preds = np.round(predictions, decimals=2)
    result_df = pd.DataFrame(rounded_preds, columns=['Probability_Class_0', 'Probability_Class_1'])
    result_df['Predicted_Class'] = np.argmax(predictions, axis=1)
    result_df['Row_Index'] = df2.index
    result_df.set_index('Row_Index', inplace=True)

    # Save probability scores
    prob_file = os.path.join(output_directory, 'probability_score.csv')
    result_df.to_csv(prob_file)
    print(f"Saved: {prob_file}")

    # Filter by threshold
    filtered_df = result_df[result_df['Probability_Class_1'] > threshold]
    if filtered_df.empty:
        print(f"No sequences found with Probability_Class_1 > {threshold}")
        return

    # Save filtered results
    filtered_file = os.path.join(output_directory, 'filtered_sequences.csv')
    filtered_df.to_csv(filtered_file)
    print(f"Saved: {filtered_file}")

    # Extract target sequences
    targets_file = os.path.join(output_directory, 'potential_targets.fasta')
    extract_sequences(modified_input, filtered_df.index, targets_file)
    print(f"Saved: {targets_file}")

    # Run DIAMOND BLAST
    print("Running DIAMOND BLAST...")
    
    # Create databases if they don't exist
    if not os.path.exists("host_protein_DB.dmnd"):
        print("Creating host protein database...")
        subprocess.run('diamond makedb --in human.fasta -d host_protein_DB', shell=True, check=True)
    
    if not os.path.exists("Essential_gene_DB.dmnd"):
        print("Creating essential gene database...")
        subprocess.run('diamond makedb --in essential.fasta -d Essential_gene_DB', shell=True, check=True)

    # BLAST (without headers for compatibility)
    host_blast = os.path.join(output_directory, "blast_against_host.tsv")
    essential_blast = os.path.join(output_directory, "blast_against_essential_genes.tsv")
    
    print("BLASTing against host proteins...")
    subprocess.run(f'diamond blastp -q {targets_file} -d host_protein_DB -o {host_blast} --very-sensitive', shell=True, check=True)
    
    print("BLASTing against essential genes...")
    subprocess.run(f'diamond blastp -q {targets_file} -d Essential_gene_DB -o {essential_blast} --very-sensitive', shell=True, check=True)

    # Add headers to BLAST results (version-independent solution)
    blast_columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 
                     'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
    header_line = '\t'.join(blast_columns) + '\n'
    
    for blast_file in [host_blast, essential_blast]:
        if os.path.exists(blast_file) and os.path.getsize(blast_file) > 0:
            with open(blast_file, 'r') as f:
                content = f.read()
            with open(blast_file, 'w') as f:
                f.write(header_line + content)
    
    print(f"✅ Analysis complete! Results saved in {output_directory}")

def extract_sequences(fasta_file, headers, output_file):
    """Extract sequences with matching headers"""
    sequences = []
    with open(fasta_file, 'r') as fasta_handle:
        for record in SeqIO.parse(fasta_handle, 'fasta'):
            if record.id in headers:
                sequences.append(record)

    with open(output_file, 'w') as output_handle:
        SeqIO.write(sequences, output_handle, 'fasta')

if __name__ == "__main__":
    main()
