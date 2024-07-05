export PATH="/hwfssz1/ST_EARTH/Reference/ST_AGRIC/APP/jdk-11.0.12+7/bin:$PATH";
export PATH="/hwfssz1/ST_EARTH/Reference/ST_AGRIC/APP/Python-3.9.6/bin/:$PATH";
for sample in `cat sample_list`
do
mkdir ${sample}
cd ${sample}
/jdfssz1/ST_EARTH/P18Z10200N0148/huangyan/08.Gene_function/Other/interproscan-5.52-86.0/interproscan.sh -goterms -dp -f tsv -cpu 2 \
 -T ./temp  \
 -i ../Protein/${sample}.fa  \
 -o  ${sample}.iprscan
cd ../
done
