cd /home/m.slagter/TONIC/fastq

print -l **/*.fastq.gz | sort > file_overview.tsv
# less file_overview.tsv

## Ensure there are no compromised samples among the listed samples
lt -ltr **/*(CF15398|CF15398)*.fastq.gz

## Amount of files
wc -l file_overview.tsv

## Amount of unique analysed samples
cat file_overview.tsv | \
  sed -E 's/fastq_files_(NORMAL|TUMOR)//g' | \
  sed -E 's/^[^_]*_[^_]*_([^_]*)_.*/\1/' | \
  sort | uniq -c | wc -l

## Total file size
du -shc **/*.fastq.gz

## Encrypt all fastqs
for f in **/*.fastq.gz; do
  echo $f
  java -jar ~/bin/EgaCryptor/EgaCryptor.jar -file $f
done

## Parallel version
# parallel -j 32 -i 'java -jar ~/bin/EgaCryptor/EgaCryptor.jar -file {}'  **/*.fastq.gz

## Upload all files
# ascp -P33001 -O33001 -QT -l300M L file.gpg  ega-box-N@fasp.ega.ebi.ac.uk:/.

for f in **/*.gpg([1]); do
  echo $f
  # ascp -P33001 -O33001 -QT -l300M L $f  m.slagter@nki.nl@fasp.ega.ebi.ac.uk:/.
done
