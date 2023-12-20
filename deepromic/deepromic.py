#!/bin/python

# File name         :deepromic.py
# description       :Classification of druggable targets using GRU model
# author            :Arittra Bhattacharjee; email: arittra.bioinfo@gmail.com 
# date of creation  :17 December 2023 to 19 December 2023 
# version           :1.0
# command           :python deepromic.py -i <path_to_input_fasta> -o <path_to_output_directory>
# required files    :human.fasta, essential.fasta and profeatx (executable file)  
# conda env create -f environment.yml -n deepromic
#**************************************************************************************************

import os
# Set the environment variable before importing TensorFlow
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"
import tensorflow as tf
# TensorFlow warnings and errors will be suppressed
tf.config.list_physical_devices("GPU")
import argparse
import subprocess
import pandas as pd
import numpy as np
import joblib
from tensorflow.keras.models import load_model
from sklearn.preprocessing import StandardScaler
from Bio import SeqIO
import warnings
warnings.filterwarnings("ignore")

print("Wellcome to DeeproMic! DeeproMic!")

# Define output directory as a global variable
output_directory = ''

def main():

    # Create ArgumentParser object
    parser = argparse.ArgumentParser(description='Description of your script.')

    # Add argument for input file
    parser.add_argument('-i', '--input', required=True, help='Input protein FASTA file (sequences with short headers are preferable)')
    
    # Add argument for probability score threshold
    parser.add_argument('-t', '--threshold', type=float, default=0.5, help=' Threshold for Probability score to filter druggable proteins (Probability_Class_1)')

    # Add argument for output directory
    parser.add_argument('-o', '--output', default='output', help='Output directory for saving generated files')

    # Parse command-line arguments
    args = parser.parse_args()

    # Access the input, threshold, and output directory from command-line arguments
    input_fasta = args.input
    threshold = args.threshold
    global output_directory
    output_directory = args.output

    # Create the output directory if it doesn't exist
    os.makedirs(output_directory, exist_ok=True)

    # Create the output directory if it doesn't exist
    os.makedirs(output_directory, exist_ok=True)

    # Modify the input FASTA file
    modified_input = input_fasta.replace('.fasta', '_modified.fasta')
    modify_command = f"awk '/^>/{{sub(/ .*/, \"\"); print; next}} {{if(NF) print}}' {input_fasta} > {modified_input}"
    subprocess.run(modify_command, shell=True, check=True)

    # Construct the shell command
    print("Converting amino acid sequences into Dipeptide Deviation from Expected Mean (DDE) features")
    command = f"./profeatx -i {modified_input} -o dde.csv -e DDE"

    # Run the shell command
    try:
        subprocess.run(command, shell=True, check=True)
        print("profeatx executed successfully.")

    except subprocess.CalledProcessError as e:
        print(f"Error: {e}")

    # Load the trained model, providing the custom loss function
    model = load_model("lstm_trained_model.h5")

    # Load the scaler model
    scaler = joblib.load("lstm_scaler_model.joblib")

    # Read the new dataset
    df = pd.read_csv("dde.csv", sep='\t')
    df2 = pd.read_csv("dde.csv", sep='\t', header=0, index_col=0)
    df = df.iloc[:, 1:-1]

    # Separate "name" feature
    X = df.values

    # Scale the features using the previously trained scaler
    X_test_scaled = scaler.transform(X)

    # Reshape the input data for GRU
    X_test_scaled_gru = np.reshape(X_test_scaled, (X_test_scaled.shape[0], 1, X_test_scaled.shape[1]))

    # Load the trained model and predict
    model = load_model("lstm_trained_model.h5")
    predictions = model.predict(X_test_scaled_gru)

    # Create a DataFrame with probability scores
    # Round the probability scores for better readability
    rounded_predictions = np.round(predictions, decimals=2)
    result_df = pd.DataFrame(rounded_predictions, columns=['Probability_Class_0', 'Probability_Class_1'])
    result_df['Predicted_Class'] = np.argmax(predictions, axis=1)  # Add predicted class based on argmax

    # Add actual labels if available
    if 'label' in df.columns:
        result_df['Actual_Class'] = df2['label'].values

    # Add row names/index from df2
    row_names = df2.index
    result_df['Row_Index'] = row_names

    # Set the 'Row_Index' column in result_df to be the same as the 'name' column in df2
    result_df.set_index('Row_Index', inplace=True)

    # Display the DataFrame
    result_df.to_csv('probability_score.csv')

    # Filter sequences based on the threshold and save to filtered_sequences.csv
    filtered_df = result_df[result_df['Probability_Class_1'] > threshold]
    filtered_df.to_csv('filtered_sequences.csv')

    # Extract sequences from the input FASTA file and save to potential_targets.fasta
    extract_sequences(input_fasta, filtered_df.index, 'potential_targets.fasta')

    # Run the diamond makedb commands
    subprocess.run('diamond makedb --in human.fasta -d host_protein_DB', shell=True, check=True)
    subprocess.run('diamond makedb --in essential.fasta -d Essential_gene_DB', shell=True, check=True)

    # Run the diamond blastp commands
    subprocess.run('diamond blastp -q potential_targets.fasta -d host_protein_DB -o blast_against_host.tsv --very-sensitive', shell=True, check=True)
    subprocess.run('diamond blastp -q potential_targets.fasta -d Essential_gene_DB -o blast_against_essential_genes.tsv --very-sensitive', shell=True, check=True)

def extract_sequences(fasta_file, headers, output_file):
    sequences = []
    with open(fasta_file, 'r') as fasta_handle:
        for record in SeqIO.parse(fasta_handle, 'fasta'):
            if record.id in headers:
                sequences.append(record)

    with open(output_file, 'w') as output_handle:
        SeqIO.write(sequences, output_handle, 'fasta')

if __name__ == "__main__":
    main()

# Move files to the output directory
files_to_move = ['dde.csv', 'probability_score.csv', 'filtered_sequences.csv', 'potential_targets.fasta', 'blast_against_host.tsv', 'blast_against_essential_genes.tsv']

for file_name in files_to_move:
    source_path = file_name
    destination_path = os.path.join(output_directory, file_name)
    
    try:
        # Move the file
        os.rename(source_path, destination_path)
        print(f"Moved {file_name} to {output_directory}")

    except FileNotFoundError:
        print(f"Error: {file_name} not found.")

    except Exception as e:
        print(f"Error: Failed to move {file_name}. Reason: {e}")
