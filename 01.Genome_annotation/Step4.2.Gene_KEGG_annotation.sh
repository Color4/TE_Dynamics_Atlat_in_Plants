scirpt_path=/jdfsbjcas1/ST_BJ/P21H28400N0232/huangyan8/Tools/KEGG/
kegg_db_path=/jdfsbjcas1/ST_BJ/P21H28400N0232/huangyan8/Tools/KEGG/
for sample in `cat sample_list`
do
/jdfsbjcas1/ST_BJ/P21H28400N0232/huangyan8/Tools/blast-2.2.26/bin/blastall  -p blastp -e 1e-05 -a 4 -m 8 -F F \
-d /jdfsbjcas1/ST_BJ/P21H28400N0232/huangyan8/ncRNA/KEGG/plant.fa \
-i ../Protein/${sample}.fa  -o ./KEGG_Blast_Result/${sample}.blast.kegg
echo ${sample} done!
perl ${scirpt_path}get_annot_info.pl -tophit 5 -topmatch 1 -id  ${scirpt_path}plant.id.annot.xls -input ./KEGG_Blast_Result/${sample}.blast.kegg -out ./KEGG_Result/${sample}.blast.kegg.xls
perl ${scirpt_path}blast2ko.pl -input /jdfsbjcas1/ST_BJ/P21H28400N0232/huangyan8/TE_Mobile/Protein/${specie}_Clean/${sample}.fa  -output  ./KEGG_Result/${sample}.ko  -blastout  ./KEGG_Blast_Result/${sample}.blast.kegg -kegg ${kegg_db_path}plant.fa #plant.fa is kegg database
perl ${scirpt_path}pathfind.pl -kegg ${kegg_db_path}plant.fa -maptitle ${scirpt_path}map_title.tab -komap ${scirpt_path}plant_ko_map.tab -fg  ./KEGG_Result/${sample}.ko -output ./KEGG_Result/${sample}.kegg.path
done
