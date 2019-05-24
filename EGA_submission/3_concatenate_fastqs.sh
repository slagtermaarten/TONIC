#!/bin/zsh
## Concatenate files that match to exactly the same biological material,
## encrypt the resulting file(s) and upload to EGA.

fastq_dir=~/TONIC/fastq
egacryptor=/home/m.slagter/bin/EgaCryptor/EgaCryptor.jar
logdir=~/TONIC/EGA_submission/logs
mkdir -p $logdir
cd $fastq_dir

## Process WXSx, 4698 is the GCF run of the RNASeq experiment
## First get all WXS files, strip to basenames using sed, strip away file specific
## information to arrive at a list of substrings that are unique to filenames that
## should be concatenated together. These IDs are called ${bn}.
ls $fastq_dir/**/*_R1_001.fastq.gz~**/4698* | sed -n -r 's\^.*/\\p' | \
  sed -n -r 's/_S[0-9]+_R1_001.fastq.gz//p' | uniq | \
  while read bn; do
    echo "Processing ${bn}"
    # bn="4698_1_CF10652_GTGGCCT"
    print -l $fastq_dir/**/${bn}_S*_R1_001.fastq.gz
    print -l $fastq_dir/**/${bn}_S*_R2_001.fastq.gz

    cat $fastq_dir/**/${bn}_S*_R1_001.fastq.gz > ${bn}_1.fastq.gz
    cat $fastq_dir/**/${bn}_S*_R2_001.fastq.gz > ${bn}_2.fastq.gz

    # ls -lhs **/${bn}*_R{1,2}_001.fastq.gz

    ## Encrypt
    java -jar $egacryptor -file ${bn}_1.fastq.gz
    java -jar $egacryptor -file ${bn}_2.fastq.gz

    ls -lhs ${bn}_{1,2}.fastq.gz

    ## Upload
    ASPERA_SCP_PASS=2ataC2BH \
      ascp -P33001 -O33001 -QT -l300M -L- \
      ${bn}_1.fastq.gz.{gpg,md5,gpg.md5} ega-box-1182@fasp.ega.ebi.ac.uk:/.

    ## Remove concatenation result files
    mv *.log $logdir
    rm -f ${bn}_{1,2}.*
  done

## Process RNAs, 4698 is the GCF run of the RNASeq experiment
ls $fastq_dir/**/4698*_R1_001.fastq.gz | sed -n -r 's\^.*/\\p' | \
  sed -n -r 's/_S[0-9]+_R1_001.fastq.gz//p' | uniq | \
  while read bn; do
    echo "Processing ${bn}"
    # bn="4698_1_CF10652_GTGGCCT"
    print -l $fastq_dir/**/${bn}_S*_R1_001.fastq.gz

    cat $fastq_dir/**/${bn}_S*_R1_001.fastq.gz > ${bn}_1.fastq.gz

    ## Encrypt
    java -jar $egacryptor -file ${bn}_1.fastq.gz

    ls -lhs ${bn}_1.fastq.gz

    ## Upload
    ASPERA_SCP_PASS=2ataC2BH \
      ascp -P33001 -O33001 -QT -l300M -L- \
      ${bn}_1.fastq.gz.{gpg,md5,gpg.md5} ega-box-1182@fasp.ega.ebi.ac.uk:/.

    ## Remove concatenation result files
    mv *.log $logdir
    rm -f ${bn}_1*
  done
