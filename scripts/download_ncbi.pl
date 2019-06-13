#!/usr/bin/perl
use LWP::Simple;
use Mozilla::CA;
getprint("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=".$ARGV[0]."&rettype=fasta");
