# _*_ coding: UTF-8 _*_

# Version information START --------------------------------------------------
VERSION_INFO = \
    """
    Author: ZHANG YUBO

    Version-01:
        2019-08-18 simulate four-taxon data set

    Version-02:
        2019-11-04 simulate four-taxon data set, with variable mutation rates, recombination rates and number of haplotype

    Version-03:
        2021-09-15 simulate four-taxon data set, with variable selection coefficients
    """
# Version information END ----------------------------------------------------

import argparse
import commands
import os
import random

###############################################################################
# read parameters
###############################################################################
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Simulate four-taxon data set")
    parser.add_argument("-o", "--Output",
                        help="Path to output directory")
    parser.add_argument("--RepNum",
                        help="Repeat all of the 1,005 evolutionary scenarios 20 times for training dataset and 2 times for test dataset",
                        default="20")
    parser.add_argument("-t", "--Time",
                        help="split time", default="1,1.5,2,3")
    parser.add_argument("-m", "--Migration",
                        help="Migration rates", default="0.05,0.1,0.5,1,5")
    parser.add_argument("-r", "--RecRate",
                        help="recombination rates", default="0.01")
    parser.add_argument("-s", "--SubRate",
                        help="substitution rates", default="0.01")
    parser.add_argument("-SC", "--SelectionCoefficient",
                        help="selection coefficients", default="0.001")
    parser.add_argument("-ST", "--SelectionTime",
                        help="selection times", default="0.5")
    parser.add_argument("-N", "--PopSize",
                        help="population size", default="1000000")
    parser.add_argument("-n", "--NumPerTaxon",
                        help="Number of haplotypes per taxon", default="8")
    parser.add_argument("-l", "--Length",
                        help="Length of sequence", default="5000")
    #parser.add_argument("--Intro",
                        #help="Types of introgression")
    parser.add_argument("--PathToMs",
                        help="Path to ms", default='/home/zhangyubo/Simulation/msdir/ms')
    parser.add_argument("--PathToSeqgen",
                        help="Path to seqgen", default='/home/zhangyubo/Simulation/Seq-Gen-1.3.4/source/seq-gen')
    parser.add_argument("--PathToMsMs",
                        help="Path to msms", default='/home/zhangyubo/Simulation/msms3.2rc-b163.jar') 
                     

    # -------------------------------------------------------------------->>>>>
    # load the paramters
    # -------------------------------------------------------------------->>>>>
    ARGS = parser.parse_args()
    OutputDir = ARGS.Output
    if not os.path.exists(OutputDir):
        os.mkdir(OutputDir)
    os.chdir(OutputDir)

    time = ARGS.Time.split(',')
    Migration = ARGS.Migration.split(',')

    RepNum = int(ARGS.RepNum) + 1
    RecRate = float(ARGS.RecRate)
    SubRate = float(ARGS.SubRate)
    PopSize = int(ARGS.PopSize)
    SelectionCoefficient = float(ARGS.SelectionCoefficient) * 2 * PopSize
    SelectionTime = ARGS.SelectionTime.split(',')

    NumPerTaxon = int(ARGS.NumPerTaxon)
    Length = int(ARGS.Length)
    #IntroType = ARGS.Intro
    PATHTOMS = ARGS.PathToMs
    PATHTOSEQGEN = ARGS.PathToSeqgen
    PATHTOMSMS = ARGS.PathToMsMs

    OutputFile = open('file_name', 'w')
    #labels = open('labels', 'w')
    TreeOutputFile = open('tree_file_name', 'w')

    Indv_com = str(5 * NumPerTaxon) + ' 1 -I 5 ' + (str(NumPerTaxon) + ' ') * 5
    if RecRate != 0:
        Rec_com = '-r %s %s' % (int(RecRate * Length), Length)
    else:
        Rec_com = ''
    Sub_com = '-l %s -s %s' % (Length, SubRate)
    '''
    intro_list = ['intro13', 'intro31', 'intro23', 'intro32', 'intro14', 'intro41', 'intro24', 'intro42', 'intro124',
                  'intro412', 'intro123', 'intro312']
    com_list = ['31', '13', '32', '23', '41', '14', '42', '24', '41', '14', '31', '13']
    index = intro_list.index(IntroType)
    '''
    str_t12 = time[0]
    str_t34 = time[1]
    str_t1234 = time[2]
    str_t12345 = time[3]

    # part 1 topo ((((1,2),3),4),5) introgression from 3 to 2
    for t_selection in SelectionTime:
        str_t_selection = str(float(t_selection) * float(str_t12))
        for str_m in Migration:
            for i in range(1, RepNum):
                file_name = str_t12 + '_' + str_t34 + '_' + str_t1234 + '_'+ str_t_selection + '_' + str_m + '_' + '32' + '_' + str(i)
                OutputFile.write('seq_' + file_name + '\n')
                TreeOutputFile.write('tree_' + file_name + '\n')
                #labels.write(str_m + '	0.0	' + str((10 - m) / 10.0) + '\n')
                if str_m == '0':
                    # topo without migration or selection
                    if SelectionCoefficient == 0.0:
                        print(
                        "java -jar %s -ms %s -ej %s 2 1 -ej %s 3 1 -ej %s 4 3 -ej %s 5 1 %s -T | tail -n +4 | grep -v // | head -n -2 > tree_%s" % (
                            PATHTOMSMS, Indv_com, str_t12, str_t1234, str_t34, str_t12345,
                            Rec_com, file_name))

                        os.system(
                        "java -jar %s -ms %s -ej %s 2 1 -ej %s 3 1 -ej %s 4 3 -ej %s 5 1 %s -T | tail -n +4 | grep -v // | head -n -2 > tree_%s" % (
                            PATHTOMSMS, Indv_com, str_t12, str_t1234, str_t34, str_t12345,
                            Rec_com, file_name))
                    # topo with selective sweep in P2
                    else:
                        print(
                        "java -jar %s -ms %s -ej %s 2 1 -ej %s 3 1 -ej %s 4 3 -ej %s 5 1 %s -Sp 0.5 -SI %s 5 0 0 0 0 0 -Smu %s -Sc 0 2 %s %s 0 -N %s -T | tail -n +4 | grep -v // | head -n -2 > tree_%s" % (
                            PATHTOMSMS, Indv_com, str_t12, str_t1234, str_t34, str_t12345,
                            Rec_com, str_t_selection, SubRate, SelectionCoefficient, SelectionCoefficient, PopSize, file_name))

                        os.system(
                        "java -jar %s -ms %s -ej %s 2 1 -ej %s 3 1 -ej %s 4 3 -ej %s 5 1 %s -Sp 0.5 -SI %s 5 0 0 0 0 0 -Smu %s -Sc 0 2 %s %s 0 -N %s -T | tail -n +4 | grep -v // | head -n -2 > tree_%s" % (
                            PATHTOMSMS, Indv_com, str_t12, str_t1234, str_t34, str_t12345,
                            Rec_com, str_t_selection, SubRate, SelectionCoefficient, SelectionCoefficient, PopSize, file_name))
                else:
                    # topo with neutral introgression from P3 to P2
                    if SelectionCoefficient == 0.0:
                        print(
                        "java -jar %s -ms %s -ej %s 2 1 -ej %s 3 1 -ej %s 4 3 -ej %s 5 1 -m 2 3 %s %s -T | tail -n +4 | grep -v // | head -n -2 > tree_%s" % (
                            PATHTOMSMS, Indv_com, str_t12, str_t1234, str_t34, str_t12345, str_m,
                            Rec_com, file_name))

                        os.system(
                        "java -jar %s -ms %s -ej %s 2 1 -ej %s 3 1 -ej %s 4 3 -ej %s 5 1 -m 2 3 %s %s -T | tail -n +4 | grep -v // | head -n -2 > tree_%s" % (
                            PATHTOMSMS, Indv_com, str_t12, str_t1234, str_t34, str_t12345, str_m,
                            Rec_com, file_name))
                    # topo with adaptive introgression from P3 to P2
                    else:
                        print(
                        "java -jar %s -ms %s -ej %s 2 1 -ej %s 3 1 -ej %s 4 3 -ej %s 5 1 -m 2 3 %s %s -SAA %s -SAa %s -Sp 0.5 -SI %s 5 0 0 1 0 0 -N %s -T | tail -n +4 | grep -v // | head -n -2 > tree_%s" % (
                            PATHTOMSMS, Indv_com, str_t12, str_t1234, str_t34, str_t12345, str_m,
                            Rec_com, SelectionCoefficient, SelectionCoefficient, str_t_selection, PopSize, file_name))

                        os.system(
                        "java -jar %s -ms %s -ej %s 2 1 -ej %s 3 1 -ej %s 4 3 -ej %s 5 1 -m 2 3 %s %s -SAA %s -SAa %s -Sp 0.5 -SI %s 5 0 0 1 0 0 -N %s -T | tail -n +4 | grep -v // | head -n -2 > tree_%s" % (
                            PATHTOMSMS, Indv_com, str_t12, str_t1234, str_t34, str_t12345, str_m,
                            Rec_com, SelectionCoefficient, SelectionCoefficient, str_t_selection, PopSize, file_name))
                partitions = commands.getoutput('wc -l tree_' + file_name)
                print("%s -q -mHKY %s -p %s < tree_%s > seq_%s" % (
                    PATHTOSEQGEN, Sub_com, partitions, file_name, file_name))

                os.system("%s -q -mHKY %s -p %s < tree_%s > seq_%s" % (
                    PATHTOSEQGEN, Sub_com, partitions, file_name, file_name))

    # part 2 topo ((((1,2),3),4),5) introgression from 2 to 3
    for t_selection in SelectionTime:
        str_t_selection = str(float(t_selection) * float(str_t12))
        for str_m in Migration:
            for i in range(1, RepNum):
                file_name = str_t12 + '_' + str_t34 + '_' + str_t1234 + '_'+ str_t_selection + '_' + str_m + '_' + '23' + '_' + str(i)
                OutputFile.write('seq_' + file_name + '\n')
                TreeOutputFile.write('tree_' + file_name + '\n')
                #labels.write(str_m + '	0.0	' + str((10 - m) / 10.0) + '\n')
                if str_m == '0':
                    if SelectionCoefficient == 0.0:   
                        os.system(
                        "java -jar %s -ms %s -ej %s 2 1 -ej %s 3 1 -ej %s 4 3 -ej %s 5 1 %s -T | tail -n +4 | grep -v // | head -n -2 > tree_%s" % (
                            PATHTOMSMS, Indv_com, str_t12, str_t1234, str_t34, str_t12345,
                            Rec_com, file_name))
                    else:
                        os.system(
                        "java -jar %s -ms %s -ej %s 2 1 -ej %s 3 1 -ej %s 4 3 -ej %s 5 1 %s -Sp 0.5 -SI %s 5 0 0 0 0 0 -Smu %s -Sc 0 3 %s %s 0 -N %s -T | tail -n +4 | grep -v // | head -n -2 > tree_%s" % (
                            PATHTOMSMS, Indv_com, str_t12, str_t1234, str_t34, str_t12345,
                        Rec_com, str_t_selection, SubRate, SelectionCoefficient, SelectionCoefficient, PopSize, file_name))
                else:
                    if SelectionCoefficient == 0.0:
                        os.system(
                        "java -jar %s -ms %s -ej %s 2 1 -ej %s 3 1 -ej %s 4 3 -ej %s 5 1 -m 3 2 %s %s -T | tail -n +4 | grep -v // | head -n -2 > tree_%s" % (
                            PATHTOMSMS, Indv_com, str_t12, str_t1234, str_t34, str_t12345, str_m,
                            Rec_com, file_name))
                    else:
                        os.system(
                        "java -jar %s -ms %s -ej %s 2 1 -ej %s 3 1 -ej %s 4 3 -ej %s 5 1 -m 3 2 %s %s -SAA %s -SAa %s -Sp 0.5 -SI %s 5 0 1 0 0 0 -N %s -T | tail -n +4 | grep -v // | head -n -2 > tree_%s" % (
                            PATHTOMSMS, Indv_com, str_t12, str_t1234, str_t34, str_t12345, str_m,
                            Rec_com, SelectionCoefficient, SelectionCoefficient, str_t_selection, PopSize, file_name))
                partitions = commands.getoutput('wc -l tree_' + file_name)
                os.system("%s -q -mHKY %s -p %s < tree_%s > seq_%s" % (
                    PATHTOSEQGEN, Sub_com, partitions, file_name, file_name))



# 2021-10-30