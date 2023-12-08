#!/bin/sh
##Data claenning
##the imported GWAS summary data should be tidies as the following colums
##the file name of micorbiome is ended with "_GCTA.assoc.raw"
##SNP     A1      A2      freq    b       se      p       N
##rs11596870      T       C       0.1413  -0.06904        0.08742 0.4298  1539
##rs868323789     C       A       0.04771 0.03995 0.1348  0.767   1530
##rs867349636     C       A       0.04756 0.03661 0.1349  0.7861  1535
##rs11593102      C       A       0.1577  0.07887 0.08044 0.327   1506
##rs1279614244    C       A       0.1577  0.07887 0.08044 0.327   1506
##rs11190359      A       G       0.07212 0.1822  0.1146  0.112   1539
##freq stands for the frequency of A1, which is A1 here.
##We utilized reference panel from 10k_reference_panel for chrom1--22 and bigcs for chrom X, we have filtered rare variant(MAF<0.001).
pheno_path=/share/home/microbiome_GCTA/
result_path=/share/home/microbiome_EAS/IV_selected_1e_5_001/
mkdir -p ${result_path}/shell
n=1
for pheno_file in `ls ${pheno_path}/*_GCTA.assoc.raw`
do
pheno_name=${pheno_file##*/}
name=${pheno_name%%_GCTA.assoc.raw}
echo ${name}
mkdir -p ${result_path}/${name}.ma
#To run the tasks faster, we split the 23 chromes into 23 separate tasks when running GCTA-COJO analysis.
for i in {1..22}
do
echo /share/home/gcta_1.94.0beta/gcta64 --bfile /share/home/10K_reference_panel_0.001.maf/10K_reference_panel_0.001.maf_add_rsid/bfile/chr${i}.CGP.beagle52.filter_MAF.0.001.addrsid --chr ${i} --cojo-file ${pheno_path}/${name}_GCTA.assoc.raw --cojo-slct --cojo-p 1e-5 --cojo-collinear 0.01 --out ${result_path}/${name}.ma/${name}chr${i} >> ${result_path}/shell/${n}_${name}_COJO.sh
done
echo /share/home/lsy_guyuqin/software/gcta_1.94.0beta/gcta64 --bfile /share/home/lsy_guyuqin/10K_reference_panel_0.001.maf/10K_reference_panel_0.001.maf_add_rsid/bfile/bigcs.X.rephasing_MAF.0.001.addrsid --chr 23 --cojo-file ${pheno_path}/${name}_GCTA.assoc.raw --cojo-slct --cojo-p 1e-5 --cojo-collinear 0.01 --out ${result_path}/${name}.ma/${name}chr23 >> ${result_path}/shell/${n}_${name}_COJO.sh
mkdir -p  ${result_path}/shell/${n}_${name}
chr=1
#The selected IVs file for every chromsome ends with ".jma.cojo".
#prepare the shell script to run in SLURM system
while read line
do
echo "#!/bin/sh
#SBATCH --nodes=1                  
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=07-00:00:00
#SBATCH --mem=5G
#SBATCH --partition=cpu
#SBATCH --job-name=${name}_${chr}
#SBATCH --output=${result_path}/shell/${n}_${name}/slurm_${name}_${chr}.sh
export PATH=/share/home/miniconda3/envs/env_gyq/bin:\$PATH
echo \"process will start at :\"
date
${line}
echo \"process end at : \"
date" > ${result_path}/shell/${n}_${name}/sbatch_${name}_${chr}.txt
let chr+=1
done < ${result_path}/shell/${n}_${name}_COJO.sh
let n+=1
done
