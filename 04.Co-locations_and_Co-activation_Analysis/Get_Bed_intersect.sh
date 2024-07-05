for specie in `cat specie_list`
do
for sample in `cat ${specie}.genomes`
do
for region in gene ncRNA TE TEpep;
do
bedtools intersect -a ./${specie}/${sample}_SV.bed   -b ./${specie}/Gene_Region/${sample}_${region}.bed  -wo  > ./Result/${specie}/${sample}.SV_vs_${region}_Region.bed;
bedtools intersect -a ./${specie}/${sample}_SV.bed   -b ./${specie}//${sample}_${region}.bed  -wo  > ./Result/${specie}/${sample}.SV_vs_${region}.bed;
done
for region in TEpep   #gene ncRNA TEpep;
do
bedtools intersect -a ./${specie}/${sample}_TE.bed   -b ./${specie}/${sample}_${region}.bed  -wo  > ./Result/${specie}/${sample}.TE_vs_${region}.bed;
done
for region in gene ncRNA;
do
bedtools intersect -a ./${specie}/${sample}_TEpep.bed  -b ./${specie}/${sample}_${region}.bed  -wo  > ./Result/${specie}/${sample}.TEpep_vs_${region}.bed;
done
done 
done



