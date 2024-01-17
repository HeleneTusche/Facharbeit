#!/bin/bash

IFS=$'\n'

#This scripts needs the Entrez Direct to work, see https://www.ncbi.nlm.nih.gov/books/NBK179288/ for download
#Otherwise, as much as possible is written in basic bash/perl in order to reduce the need for dependencies. This makes the code look a bit dirty, so feel free to adjust as desired.

#Map to all sequences longer than 4000 nt, or if no such sequences are available, map to the 5 longest ones of the taxonomy ID


mkdir joint_mapping


for l in $(cat virus_list.txt); 
do
	name=$(echo $l | cut -f 2)
	id=$(echo $l | cut -f 1)
	F="$id"_"$name"
	export F
	s="txid$id"
	echo "$F"
	cd "$F"

	a=$(head -1 "$F"_sequences.txt | cut -f 1)

	if [[ $a -gt 4000 ]] 
	then
		awk -F'\t' '$1 > 4000' "$F"_sequences.txt > "$F"_to_map.txt
	else 
		head -5 "$F"_sequences.txt > "$F"_to_map.txt
	fi

	 grep -f <(cut -f 2 "$F"_to_map.txt) "$F"_fastas.txt | while IFS=$'\t' read -r acc seq; do
        seq_end="${seq: -500}"
        l_end=${#seq_end}
        l_start=$(( ${#seq} - l_end ))
        seq_start="${seq:0:$l_start}"
        seq_end="${seq_end%%[CGT]AAAAAAAAAAAAAAAAAAAAA*}"
        echo ">$acc" >> ../joint_mapping/joint_sequences.fa
        echo "$seq_start$seq_end" >> ../joint_mapping/joint_sequences.fa
    done

    name=$(grep -w $id ../taxid_name.txt | cut -f 2)
    while IFS=$'\t' read -r acc _; do
        echo -e "$acc\t$id\t$name" 
    done < "$F"_to_map.txt >> ../acc_taxid_name.txt

    cd ..
done
