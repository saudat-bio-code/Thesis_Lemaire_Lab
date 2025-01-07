mkdir s4_mutations
cd s4_mutations
halBranchMutations /Users/saudat/Desktop/svm_all/rep_crobu/ancestor_ciona_newick_final_ros Ciona_intestinalis --snpFile  branch_cinte.bed
halBranchMutations /Users/saudat/Desktop/svm_all/rep_crobu/ancestor_ciona_newick_final_ros Ciona_robusta --snpFile  branch_crobu.bed

halBranchMutations /Users/saudat/Desktop/svm_all/rep_crobu/ancestor_ciona_newick_final_ros Ciona_robusta --refFile crobu_ins.bed --parentFile crobu_del.bed --snpFile k1.bed
halBranchMutations /Users/saudat/Desktop/svm_all/rep_crobu/ancestor_ciona_newick_final_ros Ciona_intestinalis --refFile cinte_ins.bed --parentFile cinte_del.bed --snpFile k2.bed
 
cat crobu_ins.bed branch_crobu.bed > combined.crobu
cat cinte_ins.bed branch_cinte.bed > combined.cinte
sort -k 1,1 -k2,2n combined.crobu > combined_sorted.crobu 
sort -k 1,1 -k2,2n combined.cinte > combined_sorted.cinte 

#rm k1.bed
#rm k2.bed
#make input for R
sed '1,2d' combined_sorted.cinte | sed '/^#/d' | cut -f 1-4 > combined_sorted.cinte.bed 
sed '1,2d' combined_sorted.crobu | sed '/^#/d' | cut -f 1-4 > combined_sorted.crobu.bed 



#explanation for the code
#keep in mind that del.bed has different coordinates - these are snps that are gone so you can't really show they aint there I guess.
#crobu_ins.bed: Contains insertion events that occurred specifically in the Ciona robusta lineage.
#crobu_del.bed: Contains deletion events that occurred in Ciona robusta relative to the ancestor.
#k1.bed: Contains SNPs detected in Ciona robusta. (This may overlap with branch_crobu.bed, but it
#could be filtered or formatted differently.)


#SNP files (polymorphisms) are k1.bed and branch_cinte.bed (corresponding for crobu). These flies overlap significantly but branch_cinte.bed has more info, thus I am using it for future analysis.
