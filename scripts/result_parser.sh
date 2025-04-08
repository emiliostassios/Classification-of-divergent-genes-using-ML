#!/bin/bash

#for i in $6*; do cd $i; for j in sim*; do cp $j $j.clean; done; sed -i 's/\-//g' *clean; sed -i 's/\ //g' *clean; echo $i ;python /home/etassios/ml_project/sim_rename.py; cd ../; done
#mkdir all_sims
#echo "FILENAME"
#for j in $6*; do cd $j;for i in GUT*; do awk '/>/{sub(">","&"FILENAME"_");sub(/\.seq.clean/,x)}1' $i > temp; mv  temp $i; done; cd ../;done
#for i in $6*; do cd $i; cp GUT* ../all_sims; cd ../; done 
dir=$4
###grep -A1 -h Taxon1 all_sims/GUT* > taxon1_seqs.fa; grep -A1 -h Taxon2 all_sims/GUT* > taxon2_seqs.fa; grep -A1 -h Taxon3 all_sims/GUT* > taxon3_seqs.fa; grep -A1 -h Taxon4 all_sims/GUT* > taxon4_seqs.fa
echo "TAXONS"
#find all_sims/ -name 'GUT*'|xargs grep -A1 -h Taxon1 > taxon1_seqs.fa 
#find all_sims/ -name 'GUT*'|xargs grep -A1 -h Taxon2 > taxon2_seqs.fa
#find all_sims/ -name 'GUT*'|xargs grep -A1 -h Taxon3 > taxon3_seqs.fa
#find all_sims/ -name 'GUT*'|xargs grep -A1 -h Taxon4 > taxon4_seqs.fa
#echo "MERGING"
#cat taxon*fa > simulated.fa
#sed -i 's/\-//g' simulated.fa
echo "Diamond"
diamond blastp -q $2 -d /mnt/fast_storage/sims/uhgp-50.dmnd --ultra-sensitive -o $1 -f 6 -p 18 -e 10
echo "Step 1"
awk -F"\t" '$11<0.001' $1 |awk '{print $1}' |sort |uniq > significant_genes.txt #diamond.out -> $1
echo "Step 2"
awk -F"\t" '{print $1}'  $1 |sort |uniq > temp
echo "Step 3"
###grep -vwf significant_genes.txt temp > query_genes.txt
python /mnt/fast_storage/sims/table_maker.py $1 $2 significant_genes.txt $5
echo "Step 4"
faSomeRecords $2 query_genes.txt $3
echo "Step 5"
python /mnt/fast_storage/sims/diamond_feature_extractor_mod.py $2 $5 $3 $dir/
