#!/bin/zsh
#
## Some files did not get uploaded correctly, try again with the same method
 
## Concatenate files that match to exactly the same biological material,
## encrypt the resulting file(s) and upload to EGA.

fastq_dir=~/TONIC/fastq
egacryptor=/home/m.slagter/bin/EgaCryptor/EgaCryptor.jar
logdir=~/TONIC/EGA_submission/logs
mkdir -p $logdir
cd $fastq_dir

basename_files=(
  4877_9_CF10688_CCAGTTCA
  4877_31_CF14404_TGGCTTCA 
  4877_11_CF10749_CTCAATGA
  4877_11_CF10749_CTCAATGA
  4877_7_CF10677_CGCATACA
  4877_29_CF14401_TGAAGAGA
  4877_50_CF16907_AATCCGTC
  4877_49_CF16906_AAGGACAC
  4877_48_CF16901_CAAGGAGC
  4877_48_CF16901_CAAGGAGC
  4877_45_CF15398_ACAGATTC
  4877_42_CF15187_AAGAGATC
  4877_42_CF15187_AAGAGATC
  4877_41_CF15185_AACTCACC
  4877_41_CF15185_AACTCACC
  4877_40_CF15184_ACACGACC
  4877_40_CF15184_ACACGACC
  4877_39_CF15178_CGGATTGC
  4877_37_CF15160_TTCACGCA
  4877_37_CF15160_TTCACGCA
  4877_36_CF15158_TGGTGGTA
  4877_34_CF15155_TATCAGCA
  4877_34_CF15155_TATCAGCA
  4877_33_CF15154_GAGCTGAA
  4877_33_CF15154_GAGCTGAA
  4877_32_CF15153_AATGTTGC
  4877_32_CF15153_AATGTTGC
  4877_28_CF14216_TAGGATGA
  4877_28_CF14216_TAGGATGA
  4877_21_CF11633_CGAACTTA
  4877_21_CF11633_CGAACTTA
  4877_18_CF11578_GCGAGTAA
  4877_17_CF11576_GCCACATA
  4877_17_CF11576_GCCACATA
  4877_15_CF11574_CCGTGAGA
  4877_14_CF11573_AGAGTCAA
  4877_14_CF11573_AGAGTCAA
  4877_13_CF10751_GATAGACA
  4877_13_CF10751_GATAGACA
  4877_12_CF10750_CTGAGCCA
  4877_12_CF10750_CTGAGCCA
  4877_8_CF10682_ACACAGAA
  4877_8_CF10682_ACACAGAA)

for bn in $basename_files; do
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
