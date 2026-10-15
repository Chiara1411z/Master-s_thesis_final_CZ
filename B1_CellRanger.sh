#!/bin/bash

# SINGLE-CELL CELL RANGER AUTOMATION SCRIPT

# export cellranger path
export PATH=/.../cellranger-10.0.0:$PATH

# Define Paths
# Input path: fastq.gz
input_path="/.../" 
# Outpur path: filtered_feature_bc_matrix/: matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz.; raw_feature_bc_matrix/, possorted_genome_bam.bam (+ .bai), web_summary.html, metrics_summary.csv, molecule_info.h5
output_path="/.../"
# Working path: 
SSD_path="/.../"

cd "$input_path"

# Loop through every sample directory
for d in */ ; do
  cd "$input_path"
  
  folder_name=${d%/}

  echo " Processing folder: $folder_name"

  cd "$d"

  # Extract the exact prefix Cell Ranger expects from the file itself
  thefile=$(ls *.gz | grep "_R1_" | head -n 1)
  if [ -z "$thefile" ]; then
      echo "Warning: No R1 FASTQ file found in $folder_name. Skipping..."
      continue
  fi
  sample_prefix=${thefile%%_S*}
  
  echo "--> Detected FASTQ Sample Prefix: $sample_prefix"

  # Setup Scratch Space on local fast SSD
  echo "--> Creating local SSD directory..."
  mkdir -p "$SSD_path/${folder_name}"

  echo "--> Copying FASTQ files to local SSD..."
  rsync -ahrl --progress *.gz "$SSD_path/${folder_name}/"
 

  # Run Cell Ranger Count Pipeline 
  cd "$SSD_path/${folder_name}/"
  echo "--> Launching Cell Ranger Count..."
  cellranger count --id="$sample_prefix" \
                   --transcriptome=/opt/refdata-gex-GRCh38-2024-A \
                   --fastqs="$SSD_path/${folder_name}/" \
                   --sample="$sample_prefix" \
                   --include-introns=true \
                   --localcores=24 \
                   --localmem=100 \
            	   --create-bam=true
  
  # Save Results
  mkdir -p "$output_path/${folder_name}"
  
  if [ -d "$sample_prefix" ]; then
      cd "$sample_prefix"
      rm -R SC_RNA_COUNTER_CS 2>/dev/null
      
      cp -RL "$SSD_path/${folder_name}/${sample_prefix}/"* "$output_path/${folder_name}"
  else
      echo "--> [ERROR] Cell Ranger failed to generate outputs for $sample_prefix."
  fi
  
  # Clear Local SSD Space for the next sample
  rm -R "$SSD_path/${folder_name}"
  
  cd "$input_path"
done

echo " finished"
