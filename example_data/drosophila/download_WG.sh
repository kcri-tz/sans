mkdir -p WG;
cd WG;
for f in `cat ../download_WG.list`; do wget $f; done;
for f in *.gz; do gunzip $f; done;
for f in *.fasta; do nf=$(echo $f | sed -E 's/d(...).*\.fasta/\1.fasta/g'); mv $f $nf; done;
ls *.fasta > list_WG.txt;
cd ..
