SAMPLE=XBSZ
/share/app/blast/2.11.0/bin/makeblastdb -in  ${SAMPLE}_protein.fa  -dbtype prot  -out ./DB/${SAMPLE}_db
for sample2 in ${SAMPLE} FSOY OUEC RMDB VMRD ZVWC
do
/share/app/blast/2.11.0/bin/blastp -query  ./Protein/${SAMPLE}_protein.fa   -db ./DB/${sample2}_db   -evalue 1e-10 -max_target_seqs 5  -outfmt 6  -out ${SAMPLE}_${sample2}  -num_threads 16
done
