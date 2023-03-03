# rely on twisst genomics_general dfoil vcftools IBDmix SPrime genomatnn
Scripts=/home/zhangyubo/Simulation
ERICAPath=/home/zhangyubo/ERICA

# Model Training
# four-taxon 
python ${ERICAPath}/ERICAModelTraining.py -i Training_dataset_1.txt,Training_dataset_2_2134.txt,Training_dataset_3_1324.txt,Training_dataset_4.txt,Training_dataset_5_2134.txt,Training_dataset_6_1324.txt -l Training_dataset_1_labels,Training_dataset_2_labels,Training_dataset_3_labels,Training_dataset_4_labels,Training_dataset_5_labels,Training_dataset_6_labels -p 4 -o four_taxon_model -e 3 -b 8 --Iteration 30000
# five-taxon 
python ${ERICAPath}/ERICAModelTraining.py -i Training_dataset_1.txt,Training_dataset_2.txt,Training_dataset_3.txt,Training_dataset_4.txt,Training_dataset_5.txt,Training_dataset_6.txt,Training_dataset_7.txt,Training_dataset_8.txt,Training_dataset_9.txt,Training_dataset_10.txt,Training_dataset_11.txt,Training_dataset_12.txt,Training_dataset_13.txt,Training_dataset_14.txt,Training_dataset_15.txt -l Training_dataset_1_labels,Training_dataset_2_labels,Training_dataset_3_labels,Training_dataset_4_labels,Training_dataset_5_labels,Training_dataset_6_labels,Training_dataset_7_labels,Training_dataset_8_labels,Training_dataset_9_labels,Training_dataset_10_labels,Training_dataset_11_labels,Training_dataset_12_labels,Training_dataset_13_labels,Training_dataset_14_labels,Training_dataset_15_labels -p 5 -o five_taxon_model -e 5 -b 10 --Iteration 30000

# Prediction using ERICA
# four-taxon 
python ${ERICAPath}/ERICAPrediction.py -i ${input}.txt -o ${input}_res.txt -p 4
# five-taxon
python ${ERICAPath}/ERICAPrediction.py -i ${input}.txt -o ${input}_res.txt -p 5

# Window-tree-based topology weighting
python2 ${Scripts}/genomics_general-master/raxml_sliding_windows.py -T 10 -g ${input}.geno --prefix ${input}_w500 -w 500 --windType sites --model GTRGAMMA --raxml raxmlHPC
# four-taxon
python2 ${Scripts}/twisst-master/run_twisst_parallel.py -t ${input}_w500.trees.gz -w ${input}_weights.tsv -g A 1,2,3,4,5,6,7,8 -g B 9,10,11,12,13,14,15,16 -g C 17,18,19,20,21,22,23,24 -g D 25,26,27,28,29,30,31,32 --method complete -T 10
# five-taxon
python2 ${Scripts}/twisst-master/run_twisst_parallel.py -t ${input}_w500.trees.gz -w ${input}_weights.tsv -g A 1,2,3,4,5,6,7,8 -g B 9,10,11,12,13,14,15,16 -g C 17,18,19,20,21,22,23,24 -g D 25,26,27,28,29,30,31,32 -g E 33,34,35,36,37,38,39,40 --method complete -T 10

# D and fd statistics of ABBA-BABA test
python2 ${Scripts}/genomics_general-master/ABBABABAwindows.py -g ${input}.geno -f haplo -o ${input}_${${windowsize}}.csv -w ${windowsize} -s ${windowsize} -P1 pop1 1,2,3,4,5,6,7,8 -P2 pop2 9,10,11,12,13,14,15,16 -P3 pop3 17,18,19,20,21,22,23,24 -O pop4 25,26,27,28,29,30,31,32  -T 10 --haploid 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32
python2 ${Scripts}/genomics_general-master/ABBABABAwindows.py -g ${input}.geno -f haplo -o ${input}_blocksize.csv -w ${blocksize} -s ${blocksize} -P1 pop1 1,2,3,4,5,6,7,8 -P2 pop2 9,10,11,12,13,14,15,16 -P3 pop3 17,18,19,20,21,22,23,24 -O pop4 25,26,27,28,29,30,31,32  -T 10 --haploid 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32
python2 ${Scripts}/windows_MBB.py -i ${input}_blocksize.csv -o ${input}_${windowsize}_MBB -w ${windowsize} -b ${blocksize}

# DFOIL
python ${Scripts}/sliding_window_dfoil.py -i ${input}.txt -o ${input}_dfoil.tab -w ${windowsize} -P1 1 -P2 9 -P3 17 -P4 25 -O 33

# IBDmix
# converting geno to vcf for each window
python ${Scripts}/geno2vcf.py -i ${input}.geno -o ${input}_Homo.vcf -f homo
# Gene flow from P3 to P2
vcftools --vcf ${input}_Homo.vcf --keep P2 --recode --recode-INFO-all --out ${input}_modern
python ${Scripts}/RandomSample.py -i ${input}_Homo.vcf -o ${input}_Archaic.vcf -w ${windowsize} --Indv 17,18,19,20,21,22,23,24
# Gene flow from P2 to P3
vcftools --vcf ${input}_Homo.vcf --keep P3 --recode --recode-INFO-all --out ${input}_modern
python ${Scripts}/RandomSample.py -i ${input}_Homo.vcf -o ${input}_Archaic.vcf -w ${windowsize} --Indv 9,10,11,12,13,14,15,16
${Scripts}/IBDmix/build/src/generate_gt -a ${input}_Archaic.vcf -m ${input}_modern.recode.vcf -o ${input}_IBD.gt
${Scripts}/IBDmix/build/src/ibdmix -g ${input}_IBD.gt -o ${input}_IBD.res -n Archaic -d 0 -a 0 -e 0
# error rates: 0.02
# ${Scripts}/IBDmix/build/src/ibdmix -g ${input}_IBD.gt -o ${input}_IBD.res -n Archaic -d 0 -a 0.02 -e 0.02

# SPrime
vcftools --vcf ${input}_Homo.vcf --keep P1P2 --recode --recode-INFO-all --out ${input}_P1P2 --maf 0.001
java -jar ${Scripts}/sprime-master/sprime.jar gt=${input}_P1P2.recode.vcf outgroup=P1 map=${Scripts}/AdaptiveTest.map out=${input}_sprime mu=2.5e-9 minscore=1

# genomatnn
# model A
genomatnn sim -n 10000 Nea_to_CEU.toml -j 20
genomatnn train -c Nea_to_CEU.toml
# model B
genomatnn sim -n 10000 Den_to_Papuan.toml
genomatnn train -c Den_to_Papuan.toml
