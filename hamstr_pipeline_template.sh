#!/usr/bin/env bash

#This should be a list of proteomes from which homologs will be added to the alignment
files=data/*homogID.fa

#This should be a list of alignments (can just be one gene) on which to add proteins to from the proteomes
orthogroups=/usr/local/hamstr/core_orthologs/ORTHO*

#for each alignment
for orthogroup in /usr/local/hamstr/core_orthologs/ORTHO*
do
	#This is a line I have to extract the name of the orthogroup
	name=$(echo $orthogroup | cut -c34-)
	echo $name

	#Now add a protein from each proteome if applicable
	for file in $files
	do
        	echo $file
        	#This extracts the species name so I can use it in the script later on
		species=$(echo $file | cut -c6-9)

        	#step 1
		#cd into the directory with the alignment and build a hmm for the alignment
        	cd /usr/local/hamstr/core_orthologs/$name
		hmmbuild hmm_dir/${name}.hmm ${name}.fa

        	#step 2
		#run hamstr for the alignment on the proteome. Check the hamstr help for all these options.
      		cd ~/hamstr
        	hamstr -sequence_file $file -protein -taxon $species -hmmset $name -refspec Locu_selkected -representative -central -force

        	#step 3
		#this just copies the results of hamstr to a temporary file and extracts the last 2 lines (the sequence that has been added, then adds that the the alignment and realigns using mafft, putting it back where it was so that the process can restart. It also gets rid of all the temporary files
        	cp ${orthogroup}/${name}.fa ./temp
        	tail -n 2 fa_dir*/${name}.fa >> temp
	        mafft temp > ${orthogroup}/${name}.fa
        	rm -r -f fa_dir*
	        rm -r -f hmm_search*
        	rm -r -f hamstrsearch*
		rm -f temp
	done

	#Copies the files to a handy directory
	cp ${orthogroup}/${name}.fa /Users/ab17362/OneDrive/ -/ University/ of/ Bristol/Hagfish_stuff/science/genes_from_spotted_gar_paper/hamstr_alignments 
done

