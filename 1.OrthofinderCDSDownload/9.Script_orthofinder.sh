#!/bin/bash
start=$(date +%s)
#
# do something...

## some computers require to change the ulimit...
ulimit -n 20000 
/Othofinder_program_directory/OrthoFinder/orthofinder -t 12 -f CDS_refseq_mammals_200816/Sequences_For_Orthofinder -s /Tree_for_orthofinder/treeUSED.nwk
#
end=$(date +%s)

seconds=$(echo "$end - $start" | bc)
echo $seconds' sec'

echo 'Formatted:'
awk -v t=$seconds 'BEGIN{t=int(t*1000); printf "%d:%02d:%02d\n", t/3600000, t/60000%60, t/1000%60}'
awk -v t=$seconds 'BEGIN{t=int(t*1000); printf "%d:%02d:%02d\n", t/3600000, t/60000%60, t/1000%60}' > time_Script_orthifinder.sh_1.txt
