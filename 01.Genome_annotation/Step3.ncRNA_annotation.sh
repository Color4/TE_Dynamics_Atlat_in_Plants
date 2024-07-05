for sample in `cat sample_list`
do
filename=../Genome/${sample}.fa
genome_total=988000 #kb change this to the true size of the genome
CMnumber=`grep "ACC" /jdfsbjcas1/ST_BJ/P21H28400N0232/huangyan8/Tools/Infernal/Rfam/Rfam.cm | wc -l `
Z=$(awk "BEGIN {print $genome_total*2*$CMnumber/1000}")
/jdfsbjcas1/ST_BJ/P21H28400N0232/huangyan8/Tools/Infernal/infernal-1.1.5-linux-intel-gcc/binaries/cmscan  -Z $Z --cut_ga  --rfam --nohmmonly --tblout ${sample}.tblout --fmt 2  --cpu 8 --clanin /jdfsbjcas1/ST_BJ/P21H28400N0232/huangyan8/Tools/Infernal/Rfam/Rfam.clanin /jdfsbjcas1/ST_BJ/P21H28400N0232/huangyan8/Tools/Infernal/Rfam/Rfam.cm  ${filename} > ${sample}.cmscan
awk 'BEGIN{OFS="\t";}{if(FNR==1) print "target_name\taccession\tquery_name\tquery_start\tquery_end\tstrand\tscore\tEvalue"; if(FNR>2 && $20!="=" && $0!~/^#/) print $2,$3,$4,$10,$11,$12,$17,$18; }' ${sample}.tblout > ${sample}.tblout.final.xls
done
