mkdir -p runs
dir="/data4/bio/runs-manolov/ecoli_crohn/MiraAssemblies";
for j in {1..11}
do
	printf -v i "%02d" $j 
	Rscript stat.r log/RCE${i}_Mira_RCE${i}_corrected.txt tmp/pieces/indel_pieces_RCE${i}_Mira_RCE${i}.inf pdf/RCE${i}.pdf
done

