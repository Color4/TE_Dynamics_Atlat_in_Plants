export PATH=/jdfsbjcas1/ST_BJ/P21H28400N0232/chengxin2/01.Software/Miniconda/envs/EDTA/bin:$PATH
for sample in `sample_list`
do
mkdir ${sample}
cd ${sample}
perl /jdfsbjcas1/ST_BJ/P21H28400N0232/chengxin2/01.Software/EDTA/EDTA.pl -genome ../Genome/${sample}.fa  --overwrite 0 --sensitive 0 --anno 1 --evaluate 0 -t 8 > ${sample}.V2.EDTA.run.log
if [ -e ${sample}".fa.mod.EDTA.TEanno.sum" ]
then
    rm -r ${sample}.fa.mod.EDTA.raw/*
    rm -r ${sample}.fa.mod.EDTA.final/*
    rm -r ${sample}.fa.mod.EDTA.combine/*
    rm -r ${sample}.fa.mod.EDTA.anno/*
    rm ${sample}.fa
    rm ${sample}.fa.mod
    rm ${sample}.fa.*.masked
fi
cd ../
done

