#!/bin/bash
#pslMap merged_peaks.psl genomes_fasta/Updrobu_vs_2021roscoff.psl outfile_map_upd.psl -swapMap
#pslMap merged_peaks.psl genomes_fasta/Updrobu_vs_2021roscoff.psl outfile_map_upd2.psl  -verbose=n
#axtChain -psl -linearGap=loose -faQ -faT Updrobu_vs_2021roscoff.psl HT.Ref.fasta BNKA01.1iCiinteBRoscoff_Satou.fasta Updrobu_vs_2021roscoff.chain > Updrobu_vs_2021roscoff.log
#lastal HT.Ref BNKA01.1iCiinteBRoscoff_Satou.fasta > Repeated_Updrobu_vs_2021roscoff.maf
faidx HT.Ref.fasta -i chromsizes > HT.genome.size
faidx BNKA01.1iCiinteBRoscoff_Satou.fasta -i chromsizes > BNKA01.genome.size
faidx JoinedScaffold.fasta -i chromsizes > JoinedScaffold.size

#2019HT versus joined_KH
lastal HT.Ref JoinedScaffold.fasta > Ht_ref_vs_JoinedScaffold.maf &
maf-convert psl Ht_ref_vs_JoinedScaffold.maf > Ht_ref_vs_JoinedScaffold.psl
axtChain -psl -linearGap=loose -faQ -faT Ht_ref_vs_JoinedScaffold.psl HT.Ref.fasta JoinedScaffold.fasta  Ht_ref_vs_joinedscaff.chain &
chainMergeSort Ht_ref_vs_joinedscaff.chain > Ht_vs_jsc_output_chainMergeSort
chainPreNet  Ht_vs_jsc_output_chainMergeSort HT.genome.size JoinedScaffold.size  Ht_vs_jsc_output_prenet
du -h Ht*
