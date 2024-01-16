#!/bin/bash

IFS=$'\n'

#This scripts needs the Entrez Direct to work, see https://www.ncbi.nlm.nih.gov/books/NBK179288/ for download
#Otherwise, as much as possible is written in basic bash/perl in order to reduce the need for dependencies. This makes the code look a bit dirty, so feel free to adjust as desired.

#efetch path
edirectpath="/Desktop/Facharbeit/edirect"

#Path to the annotated kaiju output file
kaiju_output="/Users/helenetusche/Desktop/Facharbeit/GSE228220_RNA_kaiju_summary_annotated.out"
#Which columns are taxonomy ID and species 
TaxidColumn="1"
SpeciesColumn="2"
SearchTerm="enterovirus"

grep -I "$SearchTerm" "$kaiju_output" |  cut -f "$TaxidColumn,$SpeciesColumn" | while IFS=$'\t' read -r id name; do
    name= "${name//[^0-9A-Za-z\n]/_}"
echo -e "$d\t$name"
done > virus_list.txt

grep -i "$SearchTerm" "$kaiju_output" | cut -f "$TaxidColumn,$SpeciesColumn" > taxid_name.txt




while IFS=$'\t' read -r id name; do
    F="$id"_"$name"
    export F
    s="txid$id"
    echo "$F"
    mkdir "$F"
    cd "$F" || exit
	$edirectpath/esearch -db nuccore -query "$s" | $edirectpath/efetch - format fasta > "$F"_fastas.fa
	
while IFS= read -r line; do
        if [[ $line =~ ^\>(.*) ]]; then
            echo -n "${BASH_REMATCH[1]}HEADERLINE"
        else
            echo -n "$line"
        fi
    done < "$F"_fastas.fa | tr -d '\n' | tr '>' '\n' | tr 'HEADERLINE' '\t' | sed 's/^\n//g' > "$F"_fastas.txt

    # Analysiere die Sequenzen und speichere die Ergebnisse in einer sortierten Textdatei
    while IFS=$'\t' read -r acc seq; do
        len="${#seq}"
        echo -e "$len\t$acc\t$seq"
    done < "$F"_fastas.txt | sort -rnk 1 > "$F"_sequences.txt

	cd ..
done
