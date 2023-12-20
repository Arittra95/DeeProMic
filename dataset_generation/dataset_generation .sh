#!/bin/bash

# activate conda environment

# Record the start time
SECONDS=0

# eval "$(conda shell.bash hook)"
conda activate bioinfo  # the environment contains cd-hit, blast, DIAMOND blast, seqkit


# Prompt the user for the FASTA file name
read -e -p "Enter the path to your pathogen FASTA file: " fasta_file

# Prompt the user for the FASTA file name
read -e -p "Enter the path to your host proteome FASTA file: " host_fasta_file

read -e -p "Enter the path to your essential gene FASTA file: " essential_gene_fasta_file

# Prompt for CD-HIT parameters
read -e -p "Enter the clustering threshold (e.g., 0.9 for 90% identity): " cluster_threshold


# Run CD-HIT with user-specified parameters
cd-hit -i "$fasta_file" -o 1_clustered_proteome.fasta -c "$cluster_threshold"


# Display a message indicating completion
echo "CD-HIT clustering completed."

seqkit stats "$fasta_file" > log1.txt
seqkit stats 1_clustered_proteome.fasta > log2.txt

sleep 3s 

echo "Starting BLASTing against host"

diamond makedb --in "$host_fasta_file" -d host_protein_DB


diamond blastp -q "$fasta_file" -d host_protein_DB -o blast_against_host.tsv --very-sensitive


awk '$11 > 1e-4 {print $1}' blast_against_host.tsv > non_homologues.txt
awk '$11 < 1e-4 {print $1}' blast_against_host.tsv > homologues.txt

awk '!seen[$0]++' non_homologues.txt > unique_non_homologues.txt
awk '!seen[$0]++' homologues.txt > unique_homologues.txt

seqkit grep -f unique_non_homologues.txt 1_clustered_proteome.fasta > 2_non_homologues_extracted_sequences.fasta 
seqkit grep -f unique_homologues.txt 1_clustered_proteome.fasta > 3_homologues_extracted_sequences.fasta 



seqkit stats 2_non_homologues_extracted_sequences.fasta  > log3.txt
seqkit stats 3_homologues_extracted_sequences.fasta  > log4.txt

echo "BLASTing against host finished"

sleep 3s 

echo "Starting BLASTing against essential gene database"


sleep 3s 

### Essential Gene Finding 

### Part 1: Non_homo_non_essential

echo "Non_homo_non_essential protein finding finding started" 


diamond makedb --in "$essential_gene_fasta_file" -d Essential_gene_DB

diamond blastp -q 2_non_homologues_extracted_sequences.fasta  -d Essential_gene_DB -o blast_non_homo.tsv --very-sensitive

awk '$11 > 1e-100 && $12 < 100 && $3 < 35 {print $1}' blast_non_homo.tsv > non_homo_non_essential.txt

awk '!seen[$0]++' non_homo_non_essential.txt > unique_non_homo_non_essential.txt

seqkit grep -f unique_non_homo_non_essential.txt 1_clustered_proteome.fasta > 3_unique_non_homo_non_essential.fasta 

seqkit stats 3_unique_non_homo_non_essential.fasta   > log4.txt


### Part 2: Non_homo_essential (Most importnat part) 

echo "Non_homo_essential protein finding finding started"


awk '$11 < 1e-100 && $12 > 100 && $3 > 35 {print $1}' blast_non_homo.tsv > non_homo_essential.txt 

awk '!seen[$0]++' non_homo_essential.txt > unique_non_homo_essential.txt

seqkit grep -f unique_non_homo_essential.txt 1_clustered_proteome.fasta > 4_unique_non_homo_essential.fasta 

seqkit stats 4_unique_non_homo_essential.fasta   > log5.txt


### Part 3: Homo_non_essential

echo "Homo_non_essential protein finding finding started"

diamond blastp -q 3_homologues_extracted_sequences.fasta  -d Essential_gene_DB -o blast_homo.tsv --very-sensitive

awk '$11 > 1e-100 && $12 < 100 && $3 < 35 {print $1}' blast_homo.tsv > homo_non_essential.txt

awk '!seen[$0]++' homo_non_essential.txt > unique_homo_non_essential.txt

seqkit grep -f unique_homo_non_essential.txt 1_clustered_proteome.fasta > 5_unique_homo_non_essential.fasta 

seqkit stats 5_unique_homo_non_essential.fasta   > log6.txt


### Part 4: Homo_essential

echo "Homo_essential protein finding finding started"

awk '$11 < 1e-100 && $12 > 100 && $3 > 35 {print $1}' blast_homo.tsv > homo_essential.txt

awk '!seen[$0]++' homo_essential.txt > unique_homo_essential.txt

seqkit grep -f unique_homo_essential.txt 1_clustered_proteome.fasta > 6_unique_homo_essential.fasta 

seqkit stats 6_unique_homo_essential.fasta    > log7.txt


echo "Essential gene sorting finished"


sleep 3s 



echo "Extracting mutual sequences"

sleep 3s 

seqkit seq --name --only-id 1_clustered_proteome.fasta > seq_ids_1_clustered_proteome.txt
seqkit seq --name --only-id 3_unique_non_homo_non_essential.fasta > seq_ids_non_homo_non_essential.txt
seqkit seq --name --only-id 4_unique_non_homo_essential.fasta > seq_ids_non_homo_essential.txt
seqkit seq --name --only-id 5_unique_homo_non_essential.fasta > seq_ids_homo_non_essential.txt
seqkit seq --name --only-id 6_unique_homo_essential.fasta > seq_ids_homo_essential.txt

sort seq_ids_non_homo_non_essential.txt seq_ids_non_homo_essential.txt seq_ids_homo_non_essential.txt seq_ids_homo_essential.txt \
| uniq -d > mutual_ids.txt


grep -v -f mutual_ids.txt seq_ids_non_homo_non_essential.txt > filtered_seq_ids_non_homo_non_essential.txt
grep -v -f mutual_ids.txt seq_ids_non_homo_essential.txt > filtered_seq_ids_non_homo_essential.txt
grep -v -f mutual_ids.txt seq_ids_homo_non_essential.txt > filtered_seq_ids_homo_non_essential.txt
grep -v -f mutual_ids.txt seq_ids_homo_essential.txt > filtered_seq_ids_homo_essential.txt 
cat mutual_ids.txt > mutual_ids_only.txt


seqkit grep -f filtered_seq_ids_non_homo_non_essential.txt 1_clustered_proteome.fasta > Final_non_homo_non_essential.fasta
seqkit grep -f filtered_seq_ids_non_homo_essential.txt 1_clustered_proteome.fasta > Final_non_homo_essential.fasta
seqkit grep -f filtered_seq_ids_homo_non_essential.txt 1_clustered_proteome.fasta > Final_homo_non_essential.fasta
seqkit grep -f filtered_seq_ids_homo_essential.txt 1_clustered_proteome.fasta > Final_homo_essential.fasta
seqkit grep -f mutual_ids_only.txt 1_clustered_proteome.fasta > Final_mutual.fasta


seqkit stats Final_non_homo_non_essential.fasta > log8.txt
seqkit stats Final_non_homo_essential.fasta > log9.txt
seqkit stats Final_homo_non_essential.fasta > log10.txt
seqkit stats Final_homo_essential.fasta > log11.txt
seqkit stats Final_mutual.fasta > log12.txt


echo "Isolating uncharacterized proteins final substraction"

# Combine sequence IDs from different files and remove duplicates
cat \
  filtered_seq_ids_non_homo_non_essential.txt \
  filtered_seq_ids_non_homo_essential.txt \
  filtered_seq_ids_homo_non_essential.txt \
  filtered_seq_ids_homo_essential.txt \
  mutual_ids_only.txt | awk '!seen[$0]++' > temp.txt

seqkit grep -v -n -f temp.txt 1_clustered_proteome.fasta > Final_uncharacterized_sequence.fasta
seqkit stats Final_uncharacterized_sequence.fasta > log13.txt



echo "Making log file"

cat log1.txt log2.txt log3.txt log4.txt log5.txt log6.txt log7.txt log8.txt log9.txt log10.txt log11.txt log12.txt log13.txt> log.txt 
rm -f log1.txt log2.txt log3.txt log4.txt log5.txt log6.txt log7.txt log8.txt log9.txt log10.txt log11.txt log12.txt log13.txt


# Record the end time
end_time=$(date +%s)


# Calculate the elapsed time
ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"

# Print the elapsed time
echo "Script execution time: $ELAPSED minutes"













