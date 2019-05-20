mkdir -p fa;
cd fa;
rm -f list.txt;
for nc in `cut -f1 ../IDs.tsv`;
do
	../../../scripts/download_ncbi.pl $nc > $nc.fasta;
	echo "$nc.fasta" >> list.txt	
done;

