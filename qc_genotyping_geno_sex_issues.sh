plink=/home/cmcuoco_sta/work/scripts_softwares/plink
INPUT=$1
ANCESTRY_LIST=$2
OUTPUT=$3

# 0. Ancestry filters
$plink --bfile $INPUT --keep $ANCESTRY_LIST --make-bed --out 0_ancestry_separation
## SNVs FILTERs
# 1. SNP: MAF 0.01, GENO 0.15, and snps-only just-acgt
$plink --bfile 0_ancestry_separation --geno 0.15 --maf 0.01 --snps-only just-acgt --make-bed --out 1_MAF_GENO
# 2. Join samples again
$plink --bfile 1_MAF_GENO --hwe 1e-6 --make-bed --out 2_HWE
# 3. Flag non-ambiguous SNP
awk '!( ($5=="A" && $6=="T") || \
        ($5=="T" && $6=="A") || \
        ($5=="G" && $6=="C") || \
        ($5=="C" && $6=="G")) {print $2}' 2_HWE.bim > 3_amb_snp_list.txt
$plink --bfile 2_HWE --extract 3_amb_snp_list.txt --make-bed --out 3_no_amb_snp
# 4. Flag duplicated SNVs (LAST SNV to use prior to prune) test
$plink --bfile 3_no_amb_snp --list-duplicate-vars suppress-first --out 4_no_dup_snps_list
# 5. Get BED, BIM, and FAM for next filters
$plink --bfile 3_no_amb_snp --exclude 4_no_dup_snps_list.dupvar --make-bed --out 5_final_SNV_filter

## SAMPLEs FILTERs ------------------------------------------------------------------------------------
# 6. Sample: mind 0.02/0.01 
$plink --bfile 5_final_SNV_filter --mind 0.01 --make-bed --out 6_mind
# 7. Prunning: 3000 1500 0.1
$plink --bfile 6_mind --indep-pairwise 3000 1500 0.1 --out 7_pruning_list
# 8. Extract SNVs from pruning
$plink --bfile 6_mind --extract 7_pruning_list.prune.in --make-bed --out 8_pruned
# 9. Het report (uses PRUNING)
$plink --bfile 8_pruned --het --out 9_het
# 10. flagging of samples with het
Rscript het_filter_3SD.R 9_het.het
# 11. Het: +- 3 * SD Removal of excluded samples (uses PRUNING)
$plink --bfile 8_pruned --keep 10_het.valid.sample --make-bed --out 11_het_filtered
# 12. Check-sex flagging (uses PRUNING)
$plink --bfile 11_het_filtered --check-sex --out 12_check_sex
# 13. flagginf samples that passed sex-check
awk '$3==0{print $1, $2}' 12_check_sex.sexcheck > 13_sex_check_pass_ambiguous.txt
grep OK 12_check_sex.sexcheck | awk '$3!=0{print $1, $2}' - > 13_sex_check_pass_with_sex.txt
cat 13_sex_check_pass_ambiguous.txt 13_sex_check_pass_with_sex.txt > 13_sex_check_pass.txt
# 14. Sample removal F > 0.8 males, F < 0.2 and XXX ambiguous (uses PRUNING)
$plink --bfile 11_het_filtered --keep 13_sex_check_pass.txt --make-bed --out 14_sex_checked
# 15. IDB flagging (uses PRUNING)
$plink --bfile 14_sex_checked --genome --out 15_relatedness
# 16. falling non-related samples
awk '$10 < 0.2 {print $1,$2}' 15_relatedness.genome | sort | uniq > 16_unrelated.txt
# 17. Related samples removal (uses PRUNING)
$plink --bfile 14_sex_checked --keep 16_unrelated.txt --make-bed --out 17_unrelated
# 18. Final files generation
$plink --bfile 5_final_SNV_filter --keep 16_unrelated.txt --make-bed --out $OUTPUT
