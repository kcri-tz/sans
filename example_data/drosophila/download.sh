mkdir -p fa;
cd fa;
for f in `cat ../download.list`; do wget $f; done;
for f in *.gz; do gunzip $f; done;
for f in *.fasta; do nf=$(echo $f | sed -E 's/d(...).*\.fasta/\1.fasta/g'); mv $f $nf; done;
ls *.fasta > list.txt;
cd ..
