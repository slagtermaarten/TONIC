cd /DATA/users/m.slagter/TONIC/EGA_submission/logs

cat *.log | grep 'Converted' | \
  sed -r 's/Converted file: //' | \
  sed -r 's/MD5 Sum://g' | \
  sed -r 's/.gz./.gz/g' | uniq | tee md5s.txt

cat *.log | grep -E 'Converted|PGP' | \
  sed -r 's/Converted file: //' | \
  sed -r 's/MD5 Sum:.+//g' | \
  sed -r 's/PGP Encrypted MD5://g' | \
  sed -r 's/.gz./.gz/g' | \
  sed -r 's/.gz./.gz/g' | \
  # sed ':a;N;$!ba;s/.gz\n/ /g' | \
  perl -p -e 's/.gz\n/.gz/' | \
  uniq | tee md5s_pgp.txt

wc -l md5s_pgp.txt
wc -l md5s.txt

# head md5s_pgp.txt
# tail md5s_pgp.txt
