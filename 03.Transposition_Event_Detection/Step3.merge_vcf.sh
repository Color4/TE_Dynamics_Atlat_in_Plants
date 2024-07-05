cd /jdfsbjcas1/ST_BJ/P21H28400N0232/huangyan8/SVMU/Step2.SVMU_to_VCF/VCF/Brassica_oleracea/Raw/
ls | while read f; 
do ls ${f}> ${f}_filelist; 
mv ${f}_filelist ${f}; 
cd ${f}; 
/jdfsbjcas1/ST_BJ/P21H28400N0232/huangyan8/Tools/SURVIVOR/Debug/SURVIVOR  merge ${f}_filelist  1000 2 1 0 0 50 ${f}.1000bp.merged.vcf
cd ../; 
done;
