for specie1 in `cat specie_list1`
do
mkdir ${specie1}
for pair in `cat ./Specie_List/${specie1}_list`
do
IFS="_" read -ra parts <<< ${pair}
specie2="${parts[1]}"
mkdir ./${specie1}/${specie1}_${specie2}
perl /jdfsbjcas1/ST_BJ/P21H28400N0232/huangyan8/Tools/DupGen_finder/DupGen_finder.pl -i  ./data -t ${specie1}   -c  ${specie2}    -o ./${specie1}/${specie1}_${specie2}
done
done
