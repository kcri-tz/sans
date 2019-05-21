mkdir -p fa;
cd fa;
rm -f list.txt;
cat ../IDs.txt | while read l;
do
	NC=$(echo $l | cut -f1 -d" "); 
	TAX=$(echo $l | cut -f2 -d " "); 
        ../../../scripts/download_ncbi.pl $NC > $TAX.fasta;
        echo "$TAX.fasta" >> list.txt
done;
cd ..;
