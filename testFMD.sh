#!/usr/bin/env bash
# testFMD.sh: test building FMD indexes with differet parameters and evaluate 
# space usage.
set -e

# Make a directory to work in
mkdir -p testFMD
cd testFMD

# Get hg38
if [ ! -d hg38 ]
then
    wget --progress=dot:giga http://hgdownload-test.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFa.tar.gz
    mkdir hg38
    tar -C hg38 -xvzf hg38.chromFa.tar.gz
    mv hg38/chroms/* hg38/
    rmdir hg38/chroms
fi

# Get HuRef by FTP globbing
if [ ! -d huRef ]
then
    mkdir huRef
    cd huRef
    
    # Make all the chromosome names. See
    # <http://stackoverflow.com/a/8789815/402891> for 0-padding.
    for CHROM in $(seq -f "%02g" 1 22) MT Un X Y
    do
        wget --progress=dot:giga ftp://ftp.ncbi.nih.gov/genomes/Homo_sapiens/CHR_${CHROM}/hs_alt_HuRef_*.fa.gz
    done
    gunzip *.gz
    cd ..
fi

# Watson is too hard to find, so download the official Korean genome
if [ ! -d koRef ]
then
    mkdir koRef
    cd koRef
    
    wget --progress=dot:giga ftp://ftp.kobic.re.kr/pub/KOBIC-KoreanGenome/fasta/*.fa.gz
    gunzip *.gz
    cd ..
fi

# Get an array of all the lists of genome files, one entry per genome.
FASTAS=("hg38/*.fa" "huRef/*.fa" "koRef/*.fa")

printf "RESULTS\tNUM_GENOMES\tBWT_BYTES\tSSA_BYTES\n"

for NUM_GENOMES in 1
do

    
    # Pull out that many of the genomes.
    SELECTED_FASTAS=${FASTAS[@]:0:${NUM_GENOMES}}
    
    echo "Genomes: ${NUM_GENOMES}"

    GENOME_FASTA="${NUM_GENOMES}genomes.fa"

    if [ ! -e ${GENOME_FASTA} ]
    then
        cat ${SELECTED_FASTAS} > ${GENOME_FASTA}
    fi

    echo "Indexing ${GENOME_FASTA}"

    time ../createIndex/createIndex ${GENOME_FASTA}-index ${GENOME_FASTA} --quiet --context 20
        
    # Check the sizes of the index.
    BWT_BYTES=$(stat -c%s ${GENOME_FASTA}-index/index.basename.bwt)
    SSA_BYTES=$(stat -c%s ${GENOME_FASTA}-index/index.basename.ssa)
    
    # Dump a grepable TSV line
    printf "RESULTS\t${NUM_GENOMES}\t${BWT_BYTES}\t${SSA_BYTES}\n"
    
    if [ -e gmon.out ]
    then
        # Make the profiling report
        gprof ../createIndex/createIndex gmon.out > ${GENOME_FASTA}.gprof.txt
    fi

done



