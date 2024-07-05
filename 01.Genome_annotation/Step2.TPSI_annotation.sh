export PERL5LIB=/jdfsbjcas1/ST_BJ/P21H28400N0232/chengxin2/01.Software/Miniconda/envs/perl/lib/perl5/site_perl/:$PERL5LIB
for sample in `cat sample_list`
do
mkdir ${sample}
cd ./${sample}
#cd ./${sample}
/jdfsbjcas1/ST_BJ/P21H28400N0232/huangyan8/Tools/TransposonPSI_08222010/transposonPSI.pl  ../Genome/${sample}.fa  nuc
rm ${sample}.fa
cd ../
done
