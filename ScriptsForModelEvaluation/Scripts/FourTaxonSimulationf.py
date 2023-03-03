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
                        help="split time", default="1,2,3")    
    parser.add_argument("-r", "--RecRate",
                        help="recombination rates", default="0.01")
    parser.add_argument("-s", "--SubRate",
                        help="substitution rates", default="0.01")
    parser.add_argument("-n", "--NumPerTaxon",
                        help="Number of haplotypes per taxon", default="8")
    parser.add_argument("-l", "--Length",
                        help="Length of sequence", default="5000")
    parser.add_argument("--PathToMs",
                        help="Path to ms", default='/home/zhangyubo/Simulation/msdir/ms')
    parser.add_argument("--PathToSeqgen",
                        help="Path to seqgen", default='/home/zhangyubo/Simulation/Seq-Gen-1.3.4/source/seq-gen')

                     

    # -------------------------------------------------------------------->>>>>
    # load the paramters
    # -------------------------------------------------------------------->>>>>
    ARGS = parser.parse_args()
    OutputDir = ARGS.Output
    if not os.path.exists(OutputDir):
        os.mkdir(OutputDir)
    os.chdir(OutputDir)

    time = ARGS.Time.split(',')

    RepNum = int(ARGS.RepNum) + 1
    RecRate = float(ARGS.RecRate)
    SubRate = float(ARGS.SubRate)
    NumPerTaxon = int(ARGS.NumPerTaxon)
    Length = int(ARGS.Length)
    PATHTOMS = ARGS.PathToMs
    PATHTOSEQGEN = ARGS.PathToSeqgen

    OutputFile = open('file_name', 'w')
    labels = open('labels', 'w')
    TreeOutputFile = open('tree_file_name', 'w')

    Indv_com = str(4 * NumPerTaxon) + ' 1 -I 4 ' + (str(NumPerTaxon) + ' ') * 4
    if RecRate != 0:
        Rec_com = '-r %s %s' % (int(RecRate * Length), Length)
    else:
        Rec_com = ''
    Sub_com = '-l %s -s %s' % (Length, SubRate)


    str_t12 = time[0]
    str_t123 = time[1]
    str_t1234 = time[2]

    # part 1 topo (((1,2),3),4) introgression from 3 to 2
    for t_intro in range(1, 10):
        for f23 in range(0, 11):
            str_t_intro = str(t_intro * float(time[0]) / 10.0)
            str_f23 = str(f23 / 10.0)
            for i in range(1, RepNum):
                file_name = str_t12 + '_' + str_t123 + '_' + str_t_intro + '_' + str_f23 + '_' + '32' + '_' + str(i)
                OutputFile.write('seq_' + file_name + '\n')
                TreeOutputFile.write('tree_' + file_name + '\n')
                labels.write(str_f23 + '    0.0    ' + str((10 - f23) / 10.0) + '\n')

                print(
                    "%s %s -ej %s 2 1 -ej %s 3 1 -ej %s 4 1 -es %s 2 %s -ej %s 5 3 %s -T | tail -n +4 | grep -v // > tree_%s" % (
                        PATHTOMS, Indv_com, str_t12, str_t123, str_t1234, str_t_intro, str_f23, str_t_intro,
                        Rec_com, file_name))
                os.system(
                    "%s %s -ej %s 2 1 -ej %s 3 1 -ej %s 4 1 -es %s 2 %s -ej %s 5 3 %s -T | tail -n +4 | grep -v // > tree_%s" % (
                        PATHTOMS, Indv_com, str_t12, str_t123, str_t1234, str_t_intro, str_f23, str_t_intro,
                        Rec_com, file_name))

                partitions = commands.getoutput('wc -l tree_' + file_name)
                #print("%s -q -mHKY %s -p %s < tree_%s > seq_%s" % (PATHTOSEQGEN, Sub_com, partitions, file_name, file_name))
                os.system("%s -q -mHKY %s -p %s < tree_%s > seq_%s" % (
                    PATHTOSEQGEN, Sub_com, partitions, file_name, file_name))

    # part 3 topo (((1,2),3),4) introgression from 2 to 3
    for t_intro in range(1, 10):
        for f23 in range(0, 11):
            str_t_intro = str(t_intro * float(time[0]) / 10.0)
            str_f23 = str(f23 / 10.0)
            for i in range(1, RepNum):
                file_name = str_t12 + '_' + str_t123 + '_' + str_t_intro + '_' + str_f23 + '_' + '23' + '_' + str(i)
                OutputFile.write('seq_' + file_name + '\n')
                TreeOutputFile.write('tree_' + file_name + '\n')
                labels.write(str_f23 + '    0.0    ' + str((10 - f23) / 10.0) + '\n')
                os.system(
                    "%s %s -ej %s 2 1 -ej %s 3 1 -ej %s 4 1 -es %s 3 %s -ej %s 5 2 %s -T | tail -n +4 | grep -v // > tree_%s" % (
                        PATHTOMS, Indv_com, str_t12, str_t123, str_t1234, str_t_intro, str_f23, str_t_intro,
                        Rec_com, file_name))

                partitions = commands.getoutput('wc -l tree_' + file_name)
                os.system("%s -q -mHKY %s -p %s < tree_%s > seq_%s" % (
                    PATHTOSEQGEN, Sub_com, partitions, file_name, file_name))



# 2021-09-22