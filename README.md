# DeeProMic (DEEp learning based druggable Protein classifier against Microorganisms)

<p align="center">
  <img src="https://github.com/Arittra95/DeeProMic/assets/57245109/932cc60c-37c0-45ab-b130-7d0710f22bd2" alt="Image description here">
</p>

A classifier to predict poteintial druggable targets of any microorganims.

## Introduction

DeeProMic is a Long Short-Term Memory (LSTM) based druggable protein classifier which has been trained on all poteintial druggable targets/ proteins (for human only) of UniProt Reference Clusters 90 (UniRef90).

It has been tested on various bacteria and eukaryotic proteins. DeeProMic takes two steps to identify druggable targets. 
First, it predicts the druggable proteins using an LSTM model. Finally, it charachterizes those proteins using Basic Local Alignment Search Tool (BLAST) against human proteome and essential proteins of Database of Essential Genes (DEG).

## How to install?

Requirements: 

1) [Miniconda](https://docs.conda.io/projects/miniconda/en/latest/) 

2) Operating system (OS): Built and tested on Ubuntu 22.04 LTS. You can try it using other Linux or MacOS based OSs. 

3) [ProFeatX](https://github.com/usubioinfo/profeatx) (if the provided profeatx dose not work in your system)

Downliad the files

```bash

cd /path/to/your/desired/directory
git clone https://github.com/Arittra95/DeeProMic.git

```
After downloading the files, you have to make a conda environment called "deepromic" using ```environment.yml``` . To do so, use these codes:

```bash
cd deepromic/
conda env create -f environment.yml -n deepromic
conda activate deepromic
unzip essential.zip

```

## How to use?


```bash
cd /path/to/your/deepromic/directory
conda activate deepromic
python deepromic.py -i <path_to_input_fasta> -o <path_to_output_directory>

```

## Options:

```
Wellcome to DeeproMic! Please provide a protein fasta file as input
usage: deepromic.py [-h] -i INPUT [-t THRESHOLD] [-o OUTPUT]

Description of your script.

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input protein FASTA file (sequences with short headers
                        are preferable)
  -t THRESHOLD, --threshold THRESHOLD
                        Threshold for Probability score to filter druggable
                        proteins (Probability_Class_1)
  -o OUTPUT, --output OUTPUT
                        Output directory for saving generated files

```

## Explanation of the Output files:

#### blast_against_essential_genes.tsv: 
BLAST outputs of ```potential_targets.fasta``` that were aligned against human.fasta. Suggestion: Avoid targets that are homologues to human proteins. 
#### blast_against_host.tsv:
BLAST outputs of ```potential_targets.fasta``` that were aligned against essential.fasta. Suggestion: Select targets that are homologues to essential proteins. 
#### dde.csv:
Dipeptide Deviation from Expected Mean features (DDE) of the ```-i INPUT``` fasta file.
#### filtered_sequences.csv:
proteins that have probability socre fpr ```Probability_Class_1``` more than the ```-t THRESHODL``` value. If you do not provide any threshold value, it will run with default value ```default=0.5```.  
#### potential_targets.fasta:
Fasta file that contains the protein sequences of ```filtered_sequences.csv```. 
#### probability_score.csv:
Probability score of each portein being druggable (Probability_Class_1) or non_druggable (Probability_Class_0).



