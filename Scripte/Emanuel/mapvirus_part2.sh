#!/bin/bash

IFS=$'\n'

#This scripts needs the Entrez Direct to work, see https://www.ncbi.nlm.nih.gov/books/NBK179288/ for download
#Otherwise, as much as possible is written in basic bash/perl in order to reduce the need for dependencies. This makes the code look a bit dirty, so feel free to adjust as desired.

#Map to all sequences longer than 4000 nt, or if no such sequences are available, map to the 5 longest ones of the taxonomy ID


mkdir joint_mapping


for l in `cat virus_list.txt`
do
	name=`echo $l | cut -f 2`
	id=`echo $l | cut -f 1`
	F=$id"_"$name
	export F
	s="txid$id"
	echo $F
	cd $F

	a=`head -1 "$F"_sequences.txt | cut -f 1`

	if [[ $a -gt 4000 ]]
	then
		cat "$F"_sequences.txt | perl -F'\t' -ane 'print if ($F[0]>4000)' > "$F"_to_map.txt
	else 
		head -5 "$F"_sequences.txt > "$F"_to_map.txt
	fi

	grep -f <(cut -f 2 "$F"_to_map.txt) "$F"_fastas.txt | perl -F'\t' -ane 'chomp($F[1]); $end=substr($F[1], -500); $l_end=length($end); $l_start=length($F[1])-$l_end; $start=substr($F[1], 0, $l_start); $end =~ s/([CGT])AAAAAAAAAAAAAAAAAAAAA*.*$/$1/g; print ">$F[0]\n$start$end\n";' >> ../joint_mapping/joint_sequences.fa

	name=`grep -w $id ../taxid_name.txt | cut -f 2`
	for acc in `cut -f 2 "$F"_to_map.txt`
	do
		echo -e "$acc\t$id\t$name" 
	done >> ../acc_taxid_name.txt
	cd ..
done
