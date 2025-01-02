#this was aquick fix I do not recommend doing it other than for this case
awk 'BEGIN {OFS="\t"};{print $1,$2,$3,$4,$3,$6}' all_conditions.bed  > all_conditions.formatted.bed 

#aligned and netted using lastz in advance
bedToPsl JoinedScaffold.size all_conditions.formatted.bed all_conditions.psl
pslMap all_conditions.psl Ht_vs_jsc_output_prenet convereted_all_conditions.psl -chainMapFile  
#gotta pick only one per all
sort -k10,10 -k2,2g convereted_all_conditions.psl | sort -u -k10,10 --merge > convereted_all_conditions_filtered.psl
#ready to go
psl2bed < convereted_all_conditions_filtered.psl > convereted_all_conditions_filtered.bed

