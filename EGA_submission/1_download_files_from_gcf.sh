#!/bin/zsh

## Maarten Slagter. 
## This script was designed to be run on a local machine which acts as a bridge
## between the HPC and the GCF's server.
## Download files one by one from the GCF, upload them to a remote location and
## remove the local copy

link="http://forge/userdata/sLpLDXtXQNrqruxTWwnV2GSuNVSSGAdgP60YoQzR/TONIC_stage1_WES/fastq_files/fastq_name_collision/"
remoteloc="m.slagter@coley:TONIC/fastq/dna/fastq_files_TUMOR"

if [[ -f index.html ]]; then
  rm index.html
fi

### Download index
wget $link

files=$(cat index.html | grep -o '<a href=['"'"'"][^"'"'"']*['"'"'"]' | sed -e 's/^<a href=["'"'"']//' -e 's/["'"'"']$//' | grep -v "mailto" | grep -v "forge" | grep -v "O=A" | grep -v "javascript")
files_array=("${(f)files}")

print "${#files} to transfer"

for f in $files_array; do
  print "Starting file: $f\n"; 
  wget "${link}${f}" && eval "rsync -avz $f --remove-source-files $remoteloc"
done

rm index.html
