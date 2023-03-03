# rely on trimmomatic bowtie2 picard-tools GATK genomics_general G-PhoCS ms Seq-Gen twisst vcf2phylip raxml tophat2
Scripts=/home/zhangyubo/Simulation
ERICAPath=/home/zhangyubo/ERICA

# Heliconius analyses
## SNP calling
### aligning re-sequencing data to genome
for i in `cat samplename`
do
java -jar ${Scripts}/trimmomatic-0.38.jar PE -threads 10 -phred33 ${i}_1.fastq ${i}_2.fastq output_paired_${i}_R1.fq output_unpaired_${i}_R1.fq output_paired_${i}_R2.fq output_unpaired_${i}_R2.fq LEADING:3 TRAILING:3 SLIDINGWINDOW:90:10 MINLEN:36 >> ${i}.log 2>&1
bowtie2 --very-sensitive-local -p 10 -x Hmel2.5 -1 output_paired_${i}_R1.fq -2 output_paired_${i}_R2.fq -S ${i}_best.sam >> ${i}.log 2>&1
samtools view -bST Hmel2.5.fa -o ${i}_noRG.bam ${i}_best.sam
java -jar ${Scripts}/picard-tools-1.96/AddOrReplaceReadGroups.jar INPUT=${i}_noRG.bam OUTPUT=${i}_std.bam SORT_ORDER=coordinate RGID=${i} RGLB=wal RGPL=illumina RGSM=${i} RGPU=none VALIDATION_STRINGENCY=LENIENT >> ${i}.log 2>&1
java -jar ${Scripts}/picard-tools-1.96/BuildBamIndex.jar INPUT=${i}_std.bam VALIDATION_STRINGENCY=LENIENT >> ${i}.log 2>&1
java -jar ${Scripts}/picard-tools-1.96/MarkDuplicates.jar INPUT=${i}_std.bam OUTPUT=${i}_std_noduplicates.bam METRICS_FILE=${i}_std.duplicate_matrics REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT >> ${i}.log 2>&1
java -jar ${Scripts}/picard-tools-1.96/BuildBamIndex.jar INPUT=${i}_std_noduplicates.bam VALIDATION_STRINGENCY=LENIENT >> ${i}.log 2>&1
java -jar ${Scripts}/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar -T RealignerTargetCreator -nt 10 -I ${i}_std_noduplicates.bam -R Hmel2.5.fa -o ${i}_forIndelAligner.intervals >> ${i}.log 2>&1
java -jar ${Scripts}/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar -T IndelRealigner -I ${i}_std_noduplicates.bam -R Hmel2.5.fa -targetIntervals ${i}_forIndelAligner.intervals -o ${i}_std_noduplicates.realign.bam -maxReads 100000 >> ${i}.log 2>&1
mv ${i}_std_noduplicates.realign.bam ${i}.bam
mv ${i}_std_noduplicates.realign.bai ${i}.bai
done

### SNP calling via GATK UnifiedGenotyper
java -jar ${Scripts}/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar -T UnifiedGenotyper -nt 12 -R Hmel2.5.fa
 -I agl_ERR260294.bam -I agl_ERS235655.bam -I agl_ERS235657.bam -I agl_SRS265431 \
 -I Hmel_ama_52.bam -I Hmel_ama_53.bam -I ERS2196343.sra.bam -I ERS2196344.sra \
 -I Htim_the_92.bam -I Htim_the_93.bam -I Htim_the_94.bam -I Htim_the_95 \
 -I Heth_nar_44.bam -I Heth_nar_45.bam -I Heth_nar_46.bam -I Heth_nar_47 \
 --heterozygosity 0.05 -stand_call_conf 50.0 -dcov 250 -o agl_ama_the_eth.vcf > agl_ama_the_eth.vcf.log 2>&1

## ERICA analyses
### (H. m. aglaope, H. m. amaryllis, H. t. thelxinoe, H. ethilla)
python ${ERICAPath}/vcf2MSA.py -i agl_ama_the_eth.vcf.gz -r Hmel2.5.fa --include Hmel201001o,Hmel202001o,Hmel203003o,Hmel204001o,Hmel205001o,Hmel206001o,Hmel207001o,Hmel208001o,Hmel209001o,Hmel210001o,Hmel211001o,Hmel212001o,Hmel213001o,Hmel214004o,Hmel215003o,Hmel216002o,Hmel217001o,Hmel218003o,Hmel219001o,Hmel220003o,Hmel221001o \
 -P1 agl_ERR260294,agl_ERS235655,agl_ERS235657,agl_SRS265431 \
 -P2 Hmel_ama_52,Hmel_ama_53,ERS2196343.sra,ERS2196344.sra \
 -P3 Htim_the_92,Htim_the_93,Htim_the_94,Htim_the_95 \
 -P4 Heth_nar_44,Heth_nar_45,Heth_nar_46,Heth_nar_47 \
 -f diplo -o agl_ama_the_eth
for i in Hmel201001o Hmel202001o Hmel203003o Hmel204001o Hmel205001o Hmel206001o Hmel207001o Hmel208001o Hmel209001o Hmel210001o Hmel211001o Hmel212001o Hmel213001o Hmel214004o Hmel215003o Hmel216002o Hmel217001o Hmel218003o Hmel219001o Hmel220003o Hmel221001o; 
do 
python ${ERICAPath}/ERICAPrediction.py -i agl_ama_the_eth_${i}.txt -o agl_ama_the_eth_${i}_res.txt -p 4
python ${ERICAPath}/ERICAVisualization.py -i agl_ama_the_eth_${i}_res.txt -o agl_ama_the_eth_${i}_50k -c ${i} -p 4 
done

## ABBA-BABA test
python2 ${Scripts}/genomics_general-master/VCF_processing/parseVCF.py -i agl_ama_the_eth.vcf.gz -o agl_ama_the_eth.geno.gz
python2 ${Scripts}/genomics_general-master/ABBABABAwindows.py -g agl_ama_the_eth.geno.gz -f phased -o agl_ama_the_eth_genome_5k_D.csv -w 5000 -s 5000 \
-P1 agl agl_ERR260294,agl_ERS235655,agl_ERS235657,agl_SRS265431 \
-P2 ama Hmel_ama_52,Hmel_ama_53,ERS2196343.sra,ERS2196344.sra \
-P3 the Htim_the_92,Htim_the_93,Htim_the_94,Htim_the_95 \
-O eth Heth_nar_44,Heth_nar_45,Heth_nar_46,Heth_nar_47
python2 ${Scripts}/windows_MBB.py -i agl_ama_the_eth_genome_5k_D.csv -o agl_ama_the_eth_genome_50k_MBB -w 50000 -b 5000

# Oryza WGA analyses
## ERICA analyses
python ${ERICAPath}/ERICAPrediction.py -i rice.conc.txt -o rice_five_pop_res.txt -p 5
python ${ERICAPath}/ERICAVisualization.py -i rice_five_pop_res.txt -o rice_five_pop_50k -p 5 

## G-PhoCS analyses
python2 ${Scripts}/run_gphocs.py -i ${Scripts}/gphocs_rice_7_pop_nomig_h.ctl -o rice_7_indv_Omeri -p A,B,C,D,E,F
python2 ${Scripts}/run_gphocs.py -i ${Scripts}/gphocs_rice_7_pop_nomig_h.ctl -o rice_7_indv_Omeri_v2 -p A,B,C,D,E,F
python2 ${Scripts}/run_gphocs.py -i ${Scripts}/gphocs_rice_7_pop_nomig_h.ctl -o rice_7_indv_Omeri_v3 -p A,B,C,D,E,F
G-PhoCS ${Scripts}/gphocs_rice_7_pop_mig_3.ctl

## simulation study of Oryza demographic model
### data simulation
mkdir DemoSimulation && cd DemoSimulation
for((i=1;i<=1000;i++)); 
do
${Scripts}/msdir/ms 5 1 -I 5 1 1 1 1 1 -n 2 1 -ej 1.25 2 1 -n 3 1 -ej 5.5 3 1 -n 4 1 -ej 5 4 3 -n 5 1 -ej 25 5 1 -en 0.025 5 17.25 -en 1.25 1 11 -en 5 3 1.5 -en 5.5 1 54 -en 25 1 17.25 -r 8 50000 -T | tail -n +4 | grep -v // > tree_${i}
echo tree_${i} >> tree_file_name
echo seq_${i} >> file_name
partitions=($(wc -l tree_${i}))
${Scripts}/Seq-Gen-1.3.4/source/seq-gen -mHKY -l 50000 -s 0.0001 -p $partitions <tree_${i} >seq_${i}
done
cd ..

### data labeling
python2 ${Scripts}/Genealogy2Weights.py -i DemoSimulation -o DemoSimulation_labels --NumOfTaxon 5  --NumPerTaxon 1

### ERICA analyses
python2 ${Scripts}/Phylip2Geno.py -i DemoSimulation -o DemoSimulation.geno --NumOfTaxon 5 -nin 1 -w 50000
python2 ${Scripts}/genomics_general-master/genoToSeq.py -g DemoSimulation.geno -s DemoSimulation.txt
sed -i '/>/d' DemoSimulation.txt
python ${ERICAPath}/ERICAPrediction.py -i DemoSimulation.txt -o DemoSimulation_res.txt -p 5
python ${ERICAPath}/ERICAVisualization.py -i DemoSimulation_res.txt -o DemoSimulation_res_50k -p 5 


# Oryza pan-genome analyses
## ERICA analyses
### (tropical japonica, temperate japonica, O. rufipogon, O. nivara, O. b arthii)
python ${ERICAPath}/vcf2MSA.py -i RicePanGenome3.vcf.gz -r Oryza_sativa.IRGSP-1.0.dna.toplevel.fa \
-P1 GP39,GP77,GP536,GP640,GP761-1,GP536,GP536,GP640 \
-P2 Nipponbare,HP14,HP44,HP48,HP314,HP103,HP45,UR28 \
-P3 W1943,W3095-2,Orufi,W3078-2,W1943,W3095-2,Orufi,W3078-2 \
-P4 W0170,W1698,W1754,W0123-1,Oniva,W1698,W0170,Oniva \
-P5 Obart -f haplo --include 1,2,3,4,5,6,7,8,9,10,11,12 -o troOsjap_tempOsjap_4Orufi_Oniva_Obart
for i in {1..12..1}; 
do 
python ${ERICAPath}/ERICAPrediction.py -i troOsjap_tempOsjap_4Orufi_Oniva_Obart_${i}.txt -o troOsjap_tempOsjap_4Orufi_Oniva_Obart_${i}_res.txt -p 5
python ${ERICAPath}/ERICAVisualization.py -i troOsjap_tempOsjap_4Orufi_Oniva_Obart_${i}_res.txt -o troOsjap_tempOsjap_4Orufi_Oniva_Obart_${i}_50k -c ${i} -p 5 
done

### (tropical japonica, temperate japonica, O. rufipogon, indica, O. b arthii)
python ${ERICAPath}/vcf2MSA.py -i RicePanGenome3.vcf.gz -r Oryza_sativa.IRGSP-1.0.dna.toplevel.fa \
-P1 GP39,GP77,GP536,GP640,GP761-1,GP536,GP536,GP640 \
-P2 Nipponbare,HP14,HP44,HP48,HP314,HP103,HP45,UR28 \
-P3 W1943,W3095-2,Orufi,W3078-2,W1943,W3095-2,Orufi,W3078-2 \
-P4 Osind,GP72,GLA4,GP22,HP119,HP274,HP362-2,HP396 \
-P5 Obart -f haplo --include 1,2,3,4,5,6,7,8,9,10,11,12 -o troOsjap_tempOsjap_4Orufi_Osind_Obart
for i in {1..12..1}; 
do 
python ${ERICAPath}/ERICAPrediction.py -i troOsjap_tempOsjap_4Orufi_Osind_Obart_${i}.txt -o troOsjap_tempOsjap_4Orufi_Osind_Obart_${i}_res.txt -p 5
python ${ERICAPath}/ERICAVisualization.py -i troOsjap_tempOsjap_4Orufi_Osind_Obart_${i}_res.txt -o troOsjap_tempOsjap_4Orufi_Osind_Obart_${i}_50k -c ${i} -p 5 
done

## genome tree
python2 vcf2phylip.py -i RicePanGenome3.vcf.gz
raxmlHPC -f a -m GTRGAMMA -p 12345 -x 12345 -# 100 -T 20 -s RicePanGenome3.phy -n RicePanGenome

## locus tree 
python2 vcf2phylip.py -i ${gene}.vcf
raxmlHPC -f a -m GTRGAMMA -p 12345 -x 12345 -# 100 -T 8 -s ${gene}.phy -n ${gene}

## Absolute divergence (dxy)
python2 ${Scripts}/genomics_general-master/VCF_processing/parseVCF.py -i RicePanGenome3.vcf.gz -o RicePanGenome3.geno.gz
python2 ${Scripts}/genomics_general-master/popgenWindows.py -g RicePanGenome3.geno.gz -o RicePanGenome_pop_50k.csv -w 50000 -s 50000 -f phased \
 -p Osjaptemp GP551,GP567,GP669,GP677,HP13-2,HP14,HP38,HP44,HP45,HP48,HP91-2,HP98,HP103,HP314,HP390,WYG7,KY131,DHX2,IL9,Koshihikari,LG31,UR28 \
 -p Osjaptro GP39,GP77,GP536,GP640,GP761-1 \
 -p Orufi Orufi,W1943,W3078-2,W3095-2,W0141,W1687,W1739,W1777,W1979,W2012 \
 -p Osind Osind,GLA4,GP3,GP22,GP51,GP72,GP772-1,HP119,HP263,HP274,HP327,HP362-2,HP383,HP396,HP407,HP486,HP492,HP517-1,HP577 \
 -p Osaus GP104,GP124,GP540,GP62,Kasalath \
 -p Oniva Oniva,W0123-1,W0170,W1698,W1754
# The absolute divergence was calculated as dxy * sites / window length

## gene expression analysis of rice populations 
for i in `cat samplename`
do
java -jar ${Scripts}/trimmomatic-0.38.jar SE -threads 16 -phred33 ${i}.fastq output_${i}.fq LEADING:3 TRAILING:3 SLIDINGWINDOW:15:30 MINLEN:36 >> trim_${i}.log 2>&1
${Scripts}/tophat-2.1.1.Linux_x86_64/tophat2 -p 20 -G MSUv7.gff3 -o ${i} MSUv7 output_${i}.fq >> ${i}_tophat.log 2>&1
done
