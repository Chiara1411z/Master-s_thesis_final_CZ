#!/bin/bash

set -euo pipefail

singularity --version

# Define Paths
# Input is Cellranger Output: possorted_genome_bam.bam (+ .bai), filtered_feature_bc_matrix/barcodes.tsv.gz
FANS_CR_Output="/.../"
# Output: clusters.tsv, cluster_genotypes.vcf, ambient_rna.txt, souporcell_minimap_tagged_sorted.bam, common_variants_covered.vcf, variants.done, consensus.npz, alt.mtx 
Souporcell_Output="/.../"

# metadata parsing:
SIF="./souporcell_release.sif"
GENOME_FA="/.../"
COMMON_VARIANTS="/.../"

THREADS=36

# Only pulls the container if it's not already sitting in the cwd
if [ ! -f "${SIF}" ]; then
    echo "Container not found locally, pulling..."
    singularity pull --arch amd64 "${SIF}" library://wheaton5/souporcell/souporcell:release
fi

# Samples + number of pooled individuals per sample (-k)
declare -A SAMPLES_K=(
  [...]=2
  [...]=3
  [...]=4

)

mkdir -p ./Input ./Output

for ID in */ do
    K=${SAMPLES_K[$ID]}

    echo "#############################################"
    echo
    echo "Working on sample ${ID} (k=${K})"
    date
    echo

    mkdir -p ./Input/${ID}
    mkdir -p ./Output/${ID}/logs

    cp -v ${FANS_CR_Output}/${ID}/outs/possorted_genome_bam.bam       ./Input/${ID}/
    cp -v ${FANS_CR_Output}/${ID}/outs/possorted_genome_bam.bam.bai   ./Input/${ID}/
    cp -v ${FANS_CR_Output}/${ID}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz ./Input/${ID}/
    gzip -d -f ./Input/${ID}/barcodes.tsv.gz

    singularity exec --cleanenv --no-home ${SIF} souporcell_pipeline.py \
        -i ./Input/${ID}/possorted_genome_bam.bam \
        -b ./Input/${ID}/barcodes.tsv \
        -f ${GENOME_FA} \
        -t ${THREADS} \
        -o ./Output/${ID} \
        -k ${K} \
        --common_variants ${COMMON_VARIANTS} \
        --skip_remap true

    mkdir -p ${Souporcell_Output}/${ID}
    mv ./Output/${ID}/* ${Souporcell_Output}/${ID}/

    rm -r ./Input/${ID}

done

echo "All samples done."
