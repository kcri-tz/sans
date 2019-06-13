mkdir -p fa;
cd fa;
rm -f list.txt;
cat ../IDs.txt | while read l;
do
	NC=$(echo $l | cut -f1 -d" "); 
	TAX=$(echo $l | cut -f2 -d " ");
	wget -O "$TAX.fasta" "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=$NC&rettype=fasta"
        echo "$TAX.fasta" >> list.txt
done;
cd ..;
