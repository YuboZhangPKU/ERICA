# rely on ms msms Seq-Gen genomics_general twisst 
Scripts=./Scripts

# four-taxon 
mkdir FourTaxon
cd FourTaxon

# 0. Data simulation of training datasets
mkdir 0Training_dataset && cd 0Training_dataset
for i in Training_dataset_1 Training_dataset_2 Training_dataset_3 Training_dataset_4 Training_dataset_5 Training_dataset_6
do
python2 ${Scripts}/FourTaxonSimulation.py -o $i --RepNum 20
done 

python2 ${Scripts}/Genealogy2Weights.py -i Training_dataset_1 -o Training_dataset_1_labels
python2 ${Scripts}/Genealogy2Weights.py -i Training_dataset_2 -o Training_dataset_2_labels --Order 2134
python2 ${Scripts}/Genealogy2Weights.py -i Training_dataset_3 -o Training_dataset_3_labels --Order 1324
python2 ${Scripts}/Genealogy2Weights.py -i Training_dataset_4 -o Training_dataset_4_labels
python2 ${Scripts}/Genealogy2Weights.py -i Training_dataset_5 -o Training_dataset_5_labels --Order 2134
python2 ${Scripts}/Genealogy2Weights.py -i Training_dataset_6 -o Training_dataset_6_labels --Order 1324
python2 ${Scripts}/ChangeSequenceOrder.py -i Training_dataset_2 -o Training_dataset_2_2134 --Order 2134
python2 ${Scripts}/ChangeSequenceOrder.py -i Training_dataset_3 -o Training_dataset_3_1324 --Order 1324
python2 ${Scripts}/ChangeSequenceOrder.py -i Training_dataset_5 -o Training_dataset_5_2134 --Order 2134
python2 ${Scripts}/ChangeSequenceOrder.py -i Training_dataset_6 -o Training_dataset_6_1324 --Order 1324

for i in Training_dataset_1 Training_dataset_2_2134 Training_dataset_3_1324 Training_dataset_4 Training_dataset_5_2134 Training_dataset_6_1324
do
python2 ${Scripts}/Phylip2Geno.py -i $i -o $i.geno
python2 ${Scripts}/genomics_general-master/genoToSeq.py -g $i.geno -s $i.txt
sed -i '/>/d' $i.txt
done
cd ..

# 1. Data simulation of test dataset D1
mkdir 1Test_dataset && cd 1Test_dataset
mkdir 0RawData && cd 0RawData
for i in Test_dataset_1 Test_dataset_2 Test_dataset_3
do
python2 ${Scripts}/FourTaxonSimulation.py -o $i --RepNum 2
done 

python2 ${Scripts}/Genealogy2Weights.py -i Test_dataset_1 -o Test_dataset_1_labels
python2 ${Scripts}/Genealogy2Weights.py -i Test_dataset_2 -o Test_dataset_2_labels --Order 2134
python2 ${Scripts}/Genealogy2Weights.py -i Test_dataset_3 -o Test_dataset_3_labels --Order 1324
python2 ${Scripts}/ChangeSequenceOrder.py -i Test_dataset_2 -o Test_dataset_2_2134 --Order 2134
python2 ${Scripts}/ChangeSequenceOrder.py -i Test_dataset_3 -o Test_dataset_3_1324 --Order 1324

for i in Test_dataset_1 Test_dataset_2_2134 Test_dataset_3_1324
do
python2 ${Scripts}/Phylip2Geno.py -i $i -o $i.geno
python2 ${Scripts}/genomics_general-master/genoToSeq.py -g $i.geno -s $i.txt
sed -i '/>/d' $i.txt
done 
cd ../..

# 2. Data simulation of test datasets D2-D7
mkdir 2VariableParameterTest && cd 2VariableParameterTest
mkdir 0RawData && cd 0RawData
# Recombination rates 
for i in 0 0.001 0.005 0.05 0.1 0.5
do
python2 ${Scripts}/FourTaxonSimulation.py -o RecRate_${i} --RepNum 2 -r ${i}
python2 ${Scripts}/Genealogy2Weights.py -i RecRate_${i} -o RecRate_${i}_labels
python2 ${Scripts}/Phylip2Geno.py -i RecRate_${i} -o RecRate_${i}.geno
python2 ${Scripts}/genomics_general-master/genoToSeq.py -g RecRate_${i}.geno -s RecRate_$i.txt
sed -i '/>/d' RecRate_$i.txt
done 

# Substitution rates
for i in 0.001 0.005 0.05 0.1 0.5
do
python2 ${Scripts}/FourTaxonSimulation.py -o SubRate_${i} --RepNum 2 -s ${i}
python2 ${Scripts}/Genealogy2Weights.py -i SubRate_${i} -o SubRate_${i}_labels
python2 ${Scripts}/Phylip2Geno.py -i SubRate_${i} -o SubRate_${i}.geno
python2 ${Scripts}/genomics_general-master/genoToSeq.py -g SubRate_${i}.geno -s SubRate_$i.txt
sed -i '/>/d' SubRate_$i.txt
done 
 
# Effective population size
for i in 0.001 0.005 0.05 0.1 0.5
do
python2 ${Scripts}/FourTaxon/FourTaxonSimulation.py -o Ne_${i} --RepNum 2 -r ${i} -s ${i}
python2 ${Scripts}/Genealogy2Weights.py -i Ne_${i} -o Ne_${i}_labels
python2 ${Scripts}/Phylip2Geno.py -i Ne_${i} -o Ne_${i}.geno
python2 ${Scripts}/genomics_general-master/genoToSeq.py -g Ne_${i}.geno -s Ne_$i.txt 
sed -i '/>/d' Ne_$i.txt
done 

# Samples per taxon
for i in 1 2 3 4 5 6 7
do
python2 ${Scripts}/FourTaxonSimulation.py -o NumPerTaxon_${i} --RepNum 2 -n ${i}
python2 ${Scripts}/Genealogy2Weights.py -i NumPerTaxon_${i} -o NumPerTaxon_${i}_labels --NumPerTaxon ${i}
python2 ${Scripts}/Phylip2Geno.py -i NumPerTaxon_${i} -o NumPerTaxon_${i}.geno -nin ${i}
python2 ${Scripts}/genomics_general-master/genoToSeq.py -g NumPerTaxon_${i}.geno -s NumPerTaxon_$i.txt
sed -i '/>/d' NumPerTaxon_$i.txt
done 

# Sequence error rates
for i in 0.001 0.01 0.02
do
python ${Scripts}/ErrorGeneration.py -i ../../1Test_dataset/0RawData/Test_dataset_1.geno -o Test_dataset_1_e${i}.geno -r ${i} -t Error
python2 ${Scripts}/genomics_general-master/genoToSeq.py -g Test_dataset_1_e${i}.geno -s Test_dataset_1_e${i}.txt
sed -i '/>/d' Test_dataset_1_e${i}.txt
done

# Sequence missing rates
for i in 0.01 0.1 0.2
do
python ${Scripts}/ErrorGeneration.py -i ../../1Test_dataset/0RawData/Test_dataset_1.geno -o Test_dataset_1_g${i}.geno -r ${i} -t Gap
python2 ${Scripts}/genomics_general-master/genoToSeq.py -g Test_dataset_1_g${i}.geno -s Test_dataset_1_g${i}.txt
sed -i '/>/d' Test_dataset_1_g${i}.txt
done
cd ../..

# 3. Data simulation of test datasets D8-D10
mkdir 3IntrogressionScenarios && cd 3IntrogressionScenarios
mkdir 0RawData && cd 0RawData
python2 ${Scripts}/FourTaxonSimulationf.py -o Proportion_123 -t 1,2,3 --RepNum 20
python2 ${Scripts}/FourTaxonSimulationf.py -o Proportion_0.5_1_1.5 -t 0.5,1,1.5 --RepNum 20
python2 ${Scripts}/FourTaxonSimulationf.py -o Proportion_0.1_0.2_0.3 -t 0.1,0.2,0.3 --RepNum 20

for i in 123 0.5_1_1.5 0.1_0.2_0.3
do
#python2 ${Scripts}/Genealogy2Weights.py -i Proportion_${i} -o Proportion_${i}_labels
python2 ${Scripts}/Phylip2Geno.py -i Proportion_${i} -o Proportion_${i}.geno
python2 ${Scripts}/genomics_general-master/genoToSeq.py -g Proportion_${i}.geno -s Proportion_${i}.txt 
sed -i '/>/d' Proportion_${i}.txt
done
cd ../..

# 4. Data simulation of test datasets D11-D19
mkdir 4AdaptiveTest && cd 4AdaptiveTest
mkdir 0RawData && cd 0RawData

mkdir AdaptiveTest_123 && cd AdaptiveTest_123
python2 ${Scripts}/FourTaxonSimulationAI.py -o AdaptiveTest_123 --RepNum 100 
python2 ${Scripts}/FourTaxonSimulationAI.py -o NeutralTest_123 --RepNum 100 --SelectionCoefficient 0 
python2 ${Scripts}/FourTaxonSimulationAI.py -o SweepTest_123 --RepNum 100 --Migration 0 
python2 ${Scripts}/FourTaxonSimulationAI.py -o NullTest_123 --RepNum 100 --Migration 0 --SelectionCoefficient 0 
cd ..

mkdir AdaptiveTest_0.5_1_1.5 && cd AdaptiveTest_0.5_1_1.5
python2 ${Scripts}/FourTaxonSimulationAI.py -o AdaptiveTest_0.5_1_1.5 --RepNum 100 -t 0.5,1,1.5 
python2 ${Scripts}/FourTaxonSimulationAI.py -o NeutralTest_0.5_1_1.5 --RepNum 100 --SelectionCoefficient 0 -t 0.5,1,1.5
python2 ${Scripts}/FourTaxonSimulationAI.py -o SweepTest_0.5_1_1.5 --RepNum 100 --Migration 0 -t 0.5,1,1.5
python2 ${Scripts}/FourTaxonSimulationAI.py -o NullTest_0.5_1_1.5 --RepNum 100 --Migration 0 --SelectionCoefficient 0 -t 0.5,1,1.5
cd ..

mkdir AdaptiveTest_0.1_0.2_0.3 && cd AdaptiveTest_0.1_0.2_0.3
python2 ${Scripts}/FourTaxonSimulationAI.py -o AdaptiveTest_0.1_0.2_0.3 --RepNum 100 -t 0.1,0.2,0.3
python2 ${Scripts}/FourTaxonSimulationAI.py -o NeutralTest_0.1_0.2_0.3 --RepNum 100 --SelectionCoefficient 0 -t 0.1,0.2,0.3
python2 ${Scripts}/FourTaxonSimulationAI.py -o SweepTest_0.1_0.2_0.3 --RepNum 100 --Migration 0 -t 0.1,0.2,0.3
python2 ${Scripts}/FourTaxonSimulationAI.py -o NullTest_0.1_0.2_0.3 --RepNum 100 --Migration 0 --SelectionCoefficient 0 -t 0.1,0.2,0.3
cd ..

mkdir AdaptiveTest_123_500k && cd AdaptiveTest_123_500k
python2 ${Scripts}/FourTaxonSimulationAI.py -o AdaptiveTest_123_500k --RepNum 100 -t 1,2,3 -N 500000 -r 0.005 -s 0.005 
python2 ${Scripts}/FourTaxonSimulationAI.py -o NeutralTest_123_500k --RepNum 100 --SelectionCoefficient 0 -t 1,2,3 -N 500000 -r 0.005 -s 0.005 
python2 ${Scripts}/FourTaxonSimulationAI.py -o SweepTest_123_500k --RepNum 100 --Migration 0 -t 1,2,3 -N 500000 -r 0.005 -s 0.005 
python2 ${Scripts}/FourTaxonSimulationAI.py -o NullTest_123_500k --RepNum 100 --Migration 0 --SelectionCoefficient 0 -t 1,2,3 -N 500000 -r 0.005 -s 0.005 
cd ..

mkdir AdaptiveTest_123_100k && cd AdaptiveTest_123_100k
python2 ${Scripts}/FourTaxonSimulationAI.py -o AdaptiveTest_123_100k --RepNum 100 -t 1,2,3 -N 100000 -r 0.001 -s 0.001 
python2 ${Scripts}/FourTaxonSimulationAI.py -o NeutralTest_123_100k --RepNum 100 --SelectionCoefficient 0 -t 1,2,3 -N 100000 -r 0.001 -s 0.001 
python2 ${Scripts}/FourTaxonSimulationAI.py -o SweepTest_123_100k --RepNum 100 --Migration 0 -t 1,2,3 -N 100000 -r 0.001 -s 0.001
python2 ${Scripts}/FourTaxonSimulationAI.py -o NullTest_123_100k --RepNum 100 --Migration 0 --SelectionCoefficient 0 -t 1,2,3 -N 100000 -r 0.001 -s 0.001 
cd ..

mkdir AdaptiveTest_123_500k_50k && cd AdaptiveTest_123_500k_50k
python2 ${Scripts}/FourTaxonSimulationAI.py -o AdaptiveTest_123_500k_50k --RepNum 100 -t 1,2,3 -N 500000 -r 0.005 -s 0.005 -l 50000 
python2 ${Scripts}/FourTaxonSimulationAI.py -o NeutralTest_123_500k_50k --RepNum 100 --SelectionCoefficient 0 -t 1,2,3 -N 500000 -r 0.005 -s 0.005 -l 50000 
python2 ${Scripts}/FourTaxonSimulationAI.py -o SweepTest_123_500k_50k --RepNum 100 --Migration 0 -t 1,2,3 -N 500000 -r 0.005 -s 0.005 -l 50000 
python2 ${Scripts}/FourTaxonSimulationAI.py -o NullTest_123_500k_50k --RepNum 100 --Migration 0 --SelectionCoefficient 0 -t 1,2,3 -N 500000 -r 0.005 -s 0.005 -l 50000 
cd ..

mkdir AdaptiveTest_123_100k_50k && cd AdaptiveTest_123_100k_50k
python2 ${Scripts}/FourTaxonSimulationAI.py -o AdaptiveTest_123_100k_50k --RepNum 100 -t 1,2,3 -N 100000 -r 0.001 -s 0.001 -l 50000
python2 ${Scripts}/FourTaxonSimulationAI.py -o NeutralTest_123_100k_50k --RepNum 100 --SelectionCoefficient 0 -t 1,2,3 -N 100000 -r 0.001 -s 0.001 -l 50000 
python2 ${Scripts}/FourTaxonSimulationAI.py -o SweepTest_123_100k_50k --RepNum 100 --Migration 0 -t 1,2,3 -N 100000 -r 0.001 -s 0.001 -l 50000
python2 ${Scripts}/FourTaxonSimulationAI.py -o NullTest_123_100k_50k --RepNum 100 --Migration 0 --SelectionCoefficient 0 -t 1,2,3 -N 100000 -r 0.001 -s 0.001 -l 50000 
cd ..


for Time in 123 0.5_1_1.5 0.1_0.2_0.3 123_500k 123_100k 123_500k_50k 123_100k_50k
do
for i in AdaptiveTest NeutralTest SweepTest NullTest
do
python2 ${Scripts}/Phylip2Geno.py -i AdaptiveTest_${Time}/${i}_${Time} -o AdaptiveTest_${Time}/${i}_${Time}.geno
python2 ${Scripts}/genomics_general-master/genoToSeq.py -g AdaptiveTest_${Time}/${i}_${Time}.geno -s AdaptiveTest_${Time}/${i}_${Time}.txt
sed -i '/>/d' AdaptiveTest_${Time}/${i}_${Time}.txt
done
done

mkdir AdaptiveTest_123_e0.02 && cd AdaptiveTest_123_e0.02
for i in AdaptiveTest NeutralTest SweepTest NullTest
do
python ${Scripts}/ErrorGeneration.py -i ../AdaptiveTest_123/${i}_123.geno -o ${i}_123_e0.02.geno -r 0.02 -t Error
python2 ${Scripts}/genomics_general-master/genoToSeq.py -g ${i}_123_e0.02.geno -s ${i}_123_e0.02.txt
sed -i '/>/d' ${i}_123_e0.02.txt
done
cd ..

mkdir AdaptiveTest_123_g0.1 && cd AdaptiveTest_123_g0.1
for i in AdaptiveTest NeutralTest SweepTest NullTest
do
python ${Scripts}/ErrorGeneration.py -i ../AdaptiveTest_123/${i}_123.geno -o ${i}_123_g0.1.geno -r 0.1 -t Gap
python2 ${Scripts}/genomics_general-master/genoToSeq.py -g ${i}_123_g0.1.geno -s ${i}_123_g0.1.txt
sed -i '/>/d' ${i}_123_g0.1.txt
done
cd ..
cd ../..
########################################################################################################################

# five-taxon
mkdir FiveTaxon
cd FiveTaxon
# 0. Data simulation of training datasets
mkdir 0Training_dataset && cd 0Training_dataset
for((i=1;i<=12;i++)); do python2 ${Scripts}/AsymmetricFiveTaxonSimulation.py -o Training_dataset_${i} --RepNum 10; done
for((i=13;i<=15;i++)); do python2 ${Scripts}/SymmetricFiveTaxonSimulation.py -o Training_dataset_${i} --RepNum 10; done

Order=('12345' '13245' '23145' '23415' '24315' '34215' '13425' '14325' '34125' '12435' '14235' '24135' '12345' '13245' '14235')
for((i=1;i<=15;i++)); do python2 ${Scripts}/Genealogy2Weights.py -i Training_dataset_${i} -o Training_dataset_${i}_labels --NumOfTaxon 5 --Order ${Order[${i}-1]}; done
for((i=2;i<=12;i++)); do python2 ${Scripts}/ChangeSequenceOrder.py -i Training_dataset_${i} -o Training_dataset_${i}_${Order[${i}-1]} --NumOfTaxon 5 --Order ${Order[${i}-1]}; done
for((i=14;i<=15;i++)); do python2 ${Scripts}/ChangeSequenceOrder.py -i Training_dataset_${i} -o Training_dataset_${i}_${Order[${i}-1]} --NumOfTaxon 5 --Order ${Order[${i}-1]}; done

for((i=1;i<=1;i++)); do python2 ${Scripts}/Phylip2Geno.py -i Training_dataset_${i} -o Training_dataset_${i}.geno --NumOfTaxon 5; done
for((i=2;i<=12;i++)); do python2 ${Scripts}/Phylip2Geno.py -i Training_dataset_${i}_${Order[${i}-1]} -o Training_dataset_${i}.geno --NumOfTaxon 5; done
for((i=13;i<=13;i++)); do python2 ${Scripts}/Phylip2Geno.py -i Training_dataset_${i} -o Training_dataset_${i}.geno --NumOfTaxon 5; done
for((i=14;i<=15;i++)); do python2 ${Scripts}/Phylip2Geno.py -i Training_dataset_${i}_${Order[${i}-1]} -o Training_dataset_${i}.geno --NumOfTaxon 5; done

for((i=1;i<=15;i++)); do  python2 ${Scripts}/genomics_general-master/genoToSeq.py -g Training_dataset_${i}.geno -s Training_dataset_${i}.txt; sed -i '/>/d' Training_dataset_${i}.txt; done
cd ..

# 1. Data simulation of test dataset D1
mkdir 1Test_dataset && cd 1Test_dataset
mkdir 0RawData && cd 0RawData
for((i=1;i<=12;i++)); do python2 ${Scripts}/AsymmetricFiveTaxonSimulation.py -o Test_dataset_${i} --RepNum 1; done
for((i=13;i<=15;i++)); do python2 ${Scripts}/SymmetricFiveTaxonSimulation.py -o Test_dataset_${i} --RepNum 1; done

for((i=1;i<=15;i++)); do python2 ${Scripts}/Genealogy2Weights.py -i Test_dataset_${i} -o Test_dataset_${i}_labels --NumOfTaxon 5 --Order ${Order[${i}-1]}; done
for((i=2;i<=12;i++)); do python2 ${Scripts}/ChangeSequenceOrder.py -i Test_dataset_${i} -o Test_dataset_${i}_${Order[${i}-1]} --NumOfTaxon 5 --Order ${Order[${i}-1]}; done
for((i=14;i<=15;i++)); do python2 ${Scripts}/ChangeSequenceOrder.py -i Test_dataset_${i} -o Test_dataset_${i}_${Order[${i}-1]} --NumOfTaxon 5 --Order ${Order[${i}-1]}; done

for((i=1;i<=1;i++)); do python2 ${Scripts}/Phylip2Geno.py -i Test_dataset_${i} -o Test_dataset_${i}.geno --NumOfTaxon 5; done
for((i=2;i<=12;i++)); do python2 ${Scripts}/Phylip2Geno.py -i Test_dataset_${i}_${Order[${i}-1]} -o Test_dataset_${i}.geno --NumOfTaxon 5; done
for((i=13;i<=13;i++)); do python2 ${Scripts}/Phylip2Geno.py -i Test_dataset_${i} -o Test_dataset_${i}.geno --NumOfTaxon 5; done
for((i=14;i<=15;i++)); do python2 ${Scripts}/Phylip2Geno.py -i Test_dataset_${i}_${Order[${i}-1]} -o Test_dataset_${i}.geno --NumOfTaxon 5; done

for((i=1;i<=15;i++)); do  python2 ${Scripts}/genomics_general-master/genoToSeq.py -g Test_dataset_${i}.geno -s Test_dataset_${i}.txt; sed -i '/>/d' Test_dataset_${i}.txt; done
cd ../..

# 2. Data simulation of test datasets D2-D7
mkdir 2VariableParameterTest && cd 2VariableParameterTest
mkdir 0RawData && cd 0RawData
# Recombination rates 
for j in 0 0.001 0.005 0.05 0.1 0.5
do
python2 ${Scripts}/AsymmetricFiveTaxonSimulation.py -o RecRate_${j}_1 -r ${j} --RepNum 1
python2 ${Scripts}/SymmetricFiveTaxonSimulation.py -o RecRate_${j}_13 -r ${j} --RepNum 1

for i in 1 13
do
python2 ${Scripts}/Genealogy2Weights.py -i RecRate_${j}_${i} -o RecRate_${j}_${i}_labels --NumOfTaxon 5
python2 ${Scripts}/Phylip2Geno.py -i RecRate_${j}_${i} -o RecRate_${j}_${i}.geno --NumOfTaxon 5
python2 ${Scripts}/genomics_general-master/genoToSeq.py -g RecRate_${j}_${i}.geno -s RecRate_${j}_${i}.txt
sed -i '/>/d' RecRate_${j}_${i}.txt
done
done

# Substitution rates
for j in 0.001 0.005 0.05 0.1 0.5
do
python2 ${Scripts}/AsymmetricFiveTaxonSimulation.py -o SubRate_${j}_1 -r ${j} --RepNum 1
python2 ${Scripts}/SymmetricFiveTaxonSimulation.py -o SubRate_${j}_13 -r ${j} --RepNum 1

for i in 1 13
do
python2 ${Scripts}/Genealogy2Weights.py -i SubRate_${j}_${i} -o SubRate_${j}_${i}_labels --NumOfTaxon 5
python2 ${Scripts}/Phylip2Geno.py -i SubRate_${j}_${i} -o SubRate_${j}_${i}.geno --NumOfTaxon 5
python2 ${Scripts}/genomics_general-master/genoToSeq.py -g SubRate_${j}_${i}.geno -s SubRate_${j}_${i}.txt
sed -i '/>/d' SubRate_${j}_${i}.txt
done
done

# Effective population size
for j in 0.001 0.005 0.05 0.1 0.5
do
python2 ${Scripts}/AsymmetricFiveTaxonSimulation.py -o Ne_${j}_1 -r ${j} -s ${j} --RepNum 1
python2 ${Scripts}/SymmetricFiveTaxonSimulation.py -o Ne_${j}_13 -r ${j} -s ${j} --RepNum 1

for i in 1 13
do
python2 ${Scripts}/Genealogy2Weights.py -i Ne_${j}_${i} -o Ne_${j}_${i}_labels --NumOfTaxon 5
python2 ${Scripts}/Phylip2Geno.py -i Ne_${j}_${i} -o Ne_${j}_${i}.geno --NumOfTaxon 5
python2 ${Scripts}/genomics_general-master/genoToSeq.py -g Ne_${j}_${i}.geno -s Ne_${j}_${i}.txt
sed -i '/>/d' Ne_${j}_${i}.txt
done
done


# Samples per taxon
for j in 1 2 3 4 5 6 7
do
python2 ${Scripts}/AsymmetricFiveTaxonSimulation.py -o NumPerTaxon_${j}_1 -n ${j} --RepNum 1
python2 ${Scripts}/SymmetricFiveTaxonSimulation.py -o NumPerTaxon_${j}_13 -n ${j} --RepNum 1

for i in 1 13
do
python2 ${Scripts}/Genealogy2Weights.py -i NumPerTaxon_${j}_${i} -o NumPerTaxon_${j}_${i}_labels --NumOfTaxon 5 --NumPerTaxon ${j}
python2 ${Scripts}/Phylip2Geno.py -i NumPerTaxon_${j}_${i} -o NumPerTaxon_${j}_${i}.geno --NumOfTaxon 5 -nin ${j}
python2 ${Scripts}/genomics_general-master/genoToSeq.py -g NumPerTaxon_${j}_${i}.geno -s NumPerTaxon_${j}_${i}.txt
sed -i '/>/d' NumPerTaxon_${j}_${i}.txt
done
done

# Sequence error rates
for i in 0.001 0.01 0.02
do
for j in 1 13
python ${Scripts}/ErrorGeneration.py -i ../../1Test_dataset/0RawData/Test_dataset_${j}.geno -o Test_dataset_${j}_e${i}.geno -r ${i} -t Error
python2 ${Scripts}/genomics_general-master/genoToSeq.py -g Test_dataset_${j}_e${i}.geno -s Test_dataset_${j}_e${i}.txt
sed -i '/>/d' Test_dataset_${j}_e${i}.txt
done
done

# Sequence missing rates
for i in 0.01 0.1 0.2
do
for j in 1 13
python ${Scripts}/ErrorGeneration.py -i ../../1Test_dataset/0RawData/Test_dataset_${j}.geno -o Test_dataset_${j}_g${i}.geno -r ${i} -t Gap
python2 ${Scripts}/genomics_general-master/genoToSeq.py -g Test_dataset_${j}_g${i}.geno -s Test_dataset_${j}_g${i}.txt
sed -i '/>/d' Test_dataset_${j}_g${i}.txt
done
done
cd ../..

# 3. Data simulation of test datasets D11-D19
mkdir 3AdaptiveTest && cd 3AdaptiveTest
mkdir 0RawData && cd 0RawData

mkdir AdaptiveTest_1_1.5_2_3 && cd AdaptiveTest_1_1.5_2_3
python2 ${Scripts}/SymmetricFiveTaxonSimulationAI.py -o AdaptiveTest_1_1.5_2_3 --RepNum 100 
python2 ${Scripts}/SymmetricFiveTaxonSimulationAI.py -o NeutralTest_1_1.5_2_3 --RepNum 100 --SelectionCoefficient 0 
python2 ${Scripts}/SymmetricFiveTaxonSimulationAI.py -o SweepTest_1_1.5_2_3 --RepNum 100 --Migration 0 
python2 ${Scripts}/SymmetricFiveTaxonSimulationAI.py -o NullTest_1_1.5_2_3 --RepNum 100 --Migration 0 --SelectionCoefficient 0 
cd ..

mkdir AdaptiveTest_0.5_0.75_1_1.5 && cd AdaptiveTest_0.5_0.75_1_1.5
python2 ${Scripts}/SymmetricFiveTaxonSimulationAI.py -o AdaptiveTest_0.5_0.75_1_1.5 --RepNum 100 -t 0.5,0.75,1,1.5 
python2 ${Scripts}/SymmetricFiveTaxonSimulationAI.py -o NeutralTest_0.5_0.75_1_1.5 --RepNum 100 --SelectionCoefficient 0 -t 0.5,0.75,1,1.5
python2 ${Scripts}/SymmetricFiveTaxonSimulationAI.py -o SweepTest_0.5_0.75_1_1.5 --RepNum 100 --Migration 0 -t 0.5,0.75,1,1.5
python2 ${Scripts}/SymmetricFiveTaxonSimulationAI.py -o NullTest_0.5_0.75_1_1.5 --RepNum 100 --Migration 0 --SelectionCoefficient 0 -t 0.5,0.75,1,1.5
cd ..

mkdir AdaptiveTest_0.1_0.15_0.2_0.3 && cd AdaptiveTest_0.1_0.15_0.2_0.3
python2 ${Scripts}/SymmetricFiveTaxonSimulationAI.py -o AdaptiveTest_0.1_0.15_0.2_0.3 --RepNum 100 -t 0.1,0.15,0.2,0.3
python2 ${Scripts}/SymmetricFiveTaxonSimulationAI.py -o NeutralTest_0.1_0.15_0.2_0.3 --RepNum 100 --SelectionCoefficient 0 -t 0.1,0.15,0.2,0.3
python2 ${Scripts}/SymmetricFiveTaxonSimulationAI.py -o SweepTest_0.1_0.15_0.2_0.3 --RepNum 100 --Migration 0 -t 0.1,0.15,0.2,0.3
python2 ${Scripts}/SymmetricFiveTaxonSimulationAI.py -o NullTest_0.1_0.15_0.2_0.3 --RepNum 100 --Migration 0 --SelectionCoefficient 0 -t 0.1,0.15,0.2,0.3
cd ..

mkdir AdaptiveTest_1_1.5_2_3_500k && cd AdaptiveTest_1_1.5_2_3_500k
python2 ${Scripts}/SymmetricFiveTaxonSimulationAI.py -o AdaptiveTest_1_1.5_2_3_500k --RepNum 100 -t 1,1.5,2,3 -N 500000 -r 0.005 -s 0.005 
python2 ${Scripts}/SymmetricFiveTaxonSimulationAI.py -o NeutralTest_1_1.5_2_3_500k --RepNum 100 --SelectionCoefficient 0 -t 1,1.5,2,3 -N 500000 -r 0.005 -s 0.005 
python2 ${Scripts}/SymmetricFiveTaxonSimulationAI.py -o SweepTest_1_1.5_2_3_500k --RepNum 100 --Migration 0 -t 1,1.5,2,3 -N 500000 -r 0.005 -s 0.005 
python2 ${Scripts}/SymmetricFiveTaxonSimulationAI.py -o NullTest_1_1.5_2_3_500k --RepNum 100 --Migration 0 --SelectionCoefficient 0 -t 1,1.5,2,3 -N 500000 -r 0.005 -s 0.005 
cd ..

mkdir AdaptiveTest_1_1.5_2_3_100k && cd AdaptiveTest_1_1.5_2_3_100k
python2 ${Scripts}/SymmetricFiveTaxonSimulationAI.py -o AdaptiveTest_1_1.5_2_3_100k --RepNum 100 -t 1,1.5,2,3 -N 100000 -r 0.001 -s 0.001 
python2 ${Scripts}/SymmetricFiveTaxonSimulationAI.py -o NeutralTest_1_1.5_2_3_100k --RepNum 100 --SelectionCoefficient 0 -t 1,1.5,2,3 -N 100000 -r 0.001 -s 0.001 
python2 ${Scripts}/SymmetricFiveTaxonSimulationAI.py -o SweepTest_1_1.5_2_3_100k --RepNum 100 --Migration 0 -t 1,1.5,2,3 -N 100000 -r 0.001 -s 0.001
python2 ${Scripts}/SymmetricFiveTaxonSimulationAI.py -o NullTest_1_1.5_2_3_100k --RepNum 100 --Migration 0 --SelectionCoefficient 0 -t 1,1.5,2,3 -N 100000 -r 0.001 -s 0.001 
cd ..

mkdir AdaptiveTest_1_1.5_2_3_500k_50k && cd AdaptiveTest_1_1.5_2_3_500k_50k
python2 ${Scripts}/SymmetricFiveTaxonSimulationAI.py -o AdaptiveTest_1_1.5_2_3_500k_50k --RepNum 100 -t 1,1.5,2,3 -N 500000 -r 0.005 -s 0.005 -l 50000 
python2 ${Scripts}/SymmetricFiveTaxonSimulationAI.py -o NeutralTest_1_1.5_2_3_500k_50k --RepNum 100 --SelectionCoefficient 0 -t 1,1.5,2,3 -N 500000 -r 0.005 -s 0.005 -l 50000 
python2 ${Scripts}/SymmetricFiveTaxonSimulationAI.py -o SweepTest_1_1.5_2_3_500k_50k --RepNum 100 --Migration 0 -t 1,1.5,2,3 -N 500000 -r 0.005 -s 0.005 -l 50000 
python2 ${Scripts}/SymmetricFiveTaxonSimulationAI.py -o NullTest_1_1.5_2_3_500k_50k --RepNum 100 --Migration 0 --SelectionCoefficient 0 -t 1,1.5,2,3 -N 500000 -r 0.005 -s 0.005 -l 50000 
cd ..

mkdir AdaptiveTest_1_1.5_2_3_100k_50k && cd AdaptiveTest_1_1.5_2_3_100k_50k
python2 ${Scripts}/SymmetricFiveTaxonSimulationAI.py -o AdaptiveTest_1_1.5_2_3_100k_50k --RepNum 100 -t 1,1.5,2,3 -N 100000 -r 0.001 -s 0.001 -l 50000
python2 ${Scripts}/SymmetricFiveTaxonSimulationAI.py -o NeutralTest_1_1.5_2_3_100k_50k --RepNum 100 --SelectionCoefficient 0 -t 1,1.5,2,3 -N 100000 -r 0.001 -s 0.001 -l 50000 
python2 ${Scripts}/SymmetricFiveTaxonSimulationAI.py -o SweepTest_1_1.5_2_3_100k_50k --RepNum 100 --Migration 0 -t 1,1.5,2,3 -N 100000 -r 0.001 -s 0.001 -l 50000
python2 ${Scripts}/SymmetricFiveTaxonSimulationAI.py -o NullTest_1_1.5_2_3_100k_50k --RepNum 100 --Migration 0 --SelectionCoefficient 0 -t 1,1.5,2,3 -N 100000 -r 0.001 -s 0.001 -l 50000 
cd ..


for Time in 1_1.5_2_3 0.5_0.75_1_1.5 0.1_0.15_0.2_0.3 1_1.5_2_3_500k 1_1.5_2_3_100k 1_1.5_2_3_500k_50k 1_1.5_2_3_100k_50k
do
for i in AdaptiveTest NeutralTest SweepTest NullTest
do
python2 ${Scripts}/Phylip2Geno.py -i AdaptiveTest_${Time}/${i}_${Time} -o AdaptiveTest_${Time}/${i}_${Time}.geno --NumOfTaxon 5
python2 ${Scripts}/genomics_general-master/genoToSeq.py -g AdaptiveTest_${Time}/${i}_${Time}.geno -s AdaptiveTest_${Time}/${i}_${Time}.txt
sed -i '/>/d' AdaptiveTest_${Time}/${i}_${Time}.txt
done
done

mkdir AdaptiveTest_1_1.5_2_3_e0.02 && cd AdaptiveTest_1_1.5_2_3_e0.02
for i in AdaptiveTest NeutralTest SweepTest NullTest
do
python ${Scripts}/ErrorGeneration.py -i ../AdaptiveTest_1_1.5_2_3/${i}_1_1.5_2_3.geno -o ${i}_1_1.5_2_3_e0.02.geno -r 0.02 -t Error
python2 ${Scripts}/genomics_general-master/genoToSeq.py -g ${i}_1_1.5_2_3_e0.02.geno -s ${i}_1_1.5_2_3_e0.02.txt
sed -i '/>/d' ${i}_1_1.5_2_3_e0.02.txt
done
cd ..

mkdir AdaptiveTest_1_1.5_2_3_g0.1 && cd AdaptiveTest_1_1.5_2_3_g0.1
for i in AdaptiveTest NeutralTest SweepTest NullTest
do
python ${Scripts}/ErrorGeneration.py -i ../AdaptiveTest_1_1.5_2_3/${i}_1_1.5_2_3.geno -o ${i}_1_1.5_2_3_g0.1.geno -r 0.1 -t Gap
python2 ${Scripts}/genomics_general-master/genoToSeq.py -g ${i}_1_1.5_2_3_g0.1.geno -s ${i}_1_1.5_2_3_g0.1.txt
sed -i '/>/d' ${i}_1_1.5_2_3_g0.1.txt
done
cd ..
cd ../..
