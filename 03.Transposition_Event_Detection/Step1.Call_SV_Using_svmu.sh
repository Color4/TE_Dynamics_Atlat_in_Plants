ref=Basmati1
for SAMPLE in `cat  Basmati1_query_sample_list`
do
/share/app/mummer/4.0.0rc1/bin/nucmer  --prefix=${SAMPLE}_${ref} ../Genomes/${ref}.fa  ../Genomes/${SAMPLE}.fa   -t 8 -c 1000
/share/app/mummer/4.0.0rc1/bin/delta-filter -q ${SAMPLE}_${ref}.delta -l 1000 -1  > ${SAMPLE}_${ref}.delta.filter
/jdfsbjcas1/ST_BJ/P21H28400N0232/huangyan8/Tools/svmu/svmu ${SAMPLE}_${ref}.delta.filter  ../Genomes/${ref}.fa   ../Genomes/${SAMPLE}.fa  h --prefix ${SAMPLE}_${ref}
done