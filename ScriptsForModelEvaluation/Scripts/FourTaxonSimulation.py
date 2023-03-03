# _*_ coding: UTF-8 _*_

# Version information START --------------------------------------------------
VERSION_INFO = \
    """
    Author: ZHANG YUBO

    Version-01:
        2019-08-18 simulate four-taxon data set

    Version-02:
        2019-11-04 simulate four-taxon data set, with variable mutation rates, recombination rates and number of haplotype

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

    # part 1 topo (((1,3),2),4)
    t13_list = [2, 4, 6, 8, 10, 12, 14, 16, 18]
    for t13 in t13_list:
        for t123 in range(t13 + 2, 22, 2):
            str_t13 = str(t13 / 10.0)
            str_t123 = str(t123 / 10.0)
            for i in range(1, RepNum):
                file_name = str_t13 + '_' + str_t123 + '_' + str(i)
                OutputFile.write('seq_' + file_name + '\n')
                TreeOutputFile.write('tree_' + file_name + '\n')
                labels.write('0\t1\t0\n')
                os.system("%s %s -ej %s 3 1 -ej %s 2 1 -ej 3 4 1 %s -T | tail -n +4 | grep -v // > tree_%s" % (
                    PATHTOMS, Indv_com, str_t13, str_t123, Rec_com, file_name))
                partitions = commands.getoutput('wc -l tree_' + file_name)
                os.system(
                    "%s -q -mHKY %s -p %s < tree_%s > seq_%s" % (PATHTOSEQGEN, Sub_com, partitions, file_name, file_name))

    # part 2 topo (((1,2),3),4) introgression from 3 to 2
    t12_list = [4, 8, 12, 16, 20]
    for t12 in t12_list:
        for t123 in range(t12 + 1, 22, 6):
            for t_intro_23 in range(1, t12, 4):
                for f23 in range(0, 10):
                    str_t12 = str(t12 / 10.0)
                    str_t123 = str(t123 / 10.0)
                    str_t_intro_23 = str(t_intro_23 / 10.0)
                    str_f23 = str(f23 / 10.0)
                    for i in range(1, RepNum):
                        file_name = str_t12 + '_' + str_t123 + '_' + str_t_intro_23 + '_' + str_f23 + '_' + '32' + '_' + str(
                            i)
                        OutputFile.write('seq_' + file_name + '\n')
                        TreeOutputFile.write('tree_' + file_name + '\n')
                        labels.write(str_f23 + '	0.0	' + str((10 - f23) / 10.0) + '\n')

                        os.system(
                            "%s %s -ej %s 2 1 -ej %s 3 1 -ej 3 4 1 -es %s 2 %s -ej %s 5 3 %s -T | tail -n +4 | grep -v // > tree_%s" % (
                                PATHTOMS, Indv_com, str_t12, str_t123, str_t_intro_23, str_f23, str_t_intro_23,
                                Rec_com, file_name))
                        partitions = commands.getoutput('wc -l tree_' + file_name)
                        os.system("%s -q -mHKY %s -p %s < tree_%s > seq_%s" % (
                            PATHTOSEQGEN, Sub_com, partitions, file_name, file_name))

    # part 3 topo (((1,2),3),4) introgression from 2 to 3
    t12_list = [4, 8, 12, 16, 20]
    for t12 in t12_list:
        for t123 in range(t12 + 1, 22, 6):
            for t_intro_23 in range(1, t12, 4):
                for f23 in range(0, 10):
                    str_t12 = str(t12 / 10.0)
                    str_t123 = str(t123 / 10.0)
                    str_t_intro_23 = str(t_intro_23 / 10.0)
                    str_f23 = str(f23 / 10.0)
                    for i in range(1, RepNum):
                        file_name = str_t12 + '_' + str_t123 + '_' + str_t_intro_23 + '_' + str_f23 + '_' + '23' + '_' + str(
                            i)
                        OutputFile.write('seq_' + file_name + '\n')
                        TreeOutputFile.write('tree_' + file_name + '\n')
                        labels.write(str_f23 + '	0.0	' + str((10 - f23) / 10.0) + '\n')

                        os.system(
                            "%s %s -ej %s 2 1 -ej %s 3 1 -ej 3 4 1 -es %s 3 %s -ej %s 5 2 %s -T | tail -n +4 | grep -v // > tree_%s" % (
                                PATHTOMS, Indv_com, str_t12, str_t123, str_t_intro_23, str_f23, str_t_intro_23,
                                Rec_com, file_name))
                        partitions = commands.getoutput('wc -l tree_' + file_name)
                        os.system("%s -q -mHKY %s -p %s < tree_%s > seq_%s" % (
                            PATHTOSEQGEN, Sub_com, partitions, file_name, file_name))

    # part 4 topo (((1,2),3),4) introgression from 3 to 1 and introgression from 3 to 2
    t12_list = [4, 8, 12, 16, 20]
    for t12 in t12_list:
        for t123 in range(t12 + 1, 22, 6):
            for t_intro_23 in range(1, t12, 4):
                for f23 in range(0, 10):
                    str_t12 = str(t12 / 10.0)
                    str_t123 = str(t123 / 10.0)
                    str_t_intro_23 = str(t_intro_23 / 10.0)
                    str_f23 = str(f23 / 10.0)
                    for i in range(1, RepNum):
                        t_intro_13 = random.randint(1, t12 - 2)
                        if t_intro_13 == t_intro_23:
                            t_intro_13 += 1
                        f13 = random.randint(0, 10)
                        str_t_intro_13 = str(t_intro_13 / 10.0)
                        str_f13 = str(f13 / 10.0)
                        f_23 = (10 - f23) * f13 / 100.0
                        f_13 = f23 * (10 - f13) / 100.0
                        if t_intro_23 < t_intro_13:
                            f_23 += (10 - f23) * (10 - f13) / 100.0
                        else:
                            f_13 += (10 - f23) * (10 - f13) / 100.0
                        f_12 = 1 - f_23 - f_13
                        file_name = str_t12 + '_' + str_t123 + '_' + str_t_intro_23 + '_' + str_t_intro_13 + '_' + str(
                            '%.2f' % f_12) + '_' + str(f_13) + '_' + str(f_23) + '_' + '3132' + '_' + str(i)
                        OutputFile.write('seq_' + file_name + '\n')
                        TreeOutputFile.write('tree_' + file_name + '\n')
                        labels.write(str('%.2f' % f_12) + '	' + str(f_13) + '	' + str(f_23) + '\n')
                        if t_intro_23 <= t_intro_13:
                            os.system(
                                "%s %s -ej %s 2 1 -ej %s 3 1 -ej 3 4 1 -es %s 2 %s -ej %s 5 3 -es %s 1 %s -ej %s 6 3 %s -T | tail -n +4 | grep -v // > tree_%s" % (
                                    PATHTOMS, Indv_com, str_t12, str_t123, str_t_intro_23, str_f23, str_t_intro_23,
                                    str_t_intro_13, str_f13, str_t_intro_13, Rec_com, file_name))
                        else:
                            os.system(
                                "%s %s -ej %s 2 1 -ej %s 3 1 -ej 3 4 1 -es %s 1 %s -ej %s 5 3 -es %s 2 %s -ej %s 6 3 %s -T | tail -n +4 | grep -v // > tree_%s" % (
                                    PATHTOMS, Indv_com, str_t12, str_t123, str_t_intro_13, str_f13, str_t_intro_13,
                                    str_t_intro_23, str_f23, str_t_intro_23, Rec_com, file_name))

                        partitions = commands.getoutput('wc -l tree_' + file_name)
                        os.system("%s -q -mHKY %s -p %s < tree_%s > seq_%s" % (
                            PATHTOSEQGEN, Sub_com, partitions, file_name, file_name))

    # part 5 topo (((1,2),3),4) introgression from 1 to 3 and introgression from 2 to 3
    t12_list = [4, 8, 12, 16, 20]
    for t12 in t12_list:
        for t123 in range(t12 + 1, 22, 6):
            for t_intro_23 in range(1, t12, 4):
                for f23 in range(0, 10):
                    str_t12 = str(t12 / 10.0)
                    str_t123 = str(t123 / 10.0)
                    str_t_intro_23 = str(t_intro_23 / 10.0)
                    str_f23 = str(f23 / 10.0)
                    for i in range(1, RepNum):
                        t_intro_13 = random.randint(1, t12 - 2)
                        if t_intro_13 == t_intro_23:
                            t_intro_13 += 1
                        f13 = random.randint(0, 10)
                        str_t_intro_13 = str(t_intro_13 / 10.0)
                        str_f13 = str(f13 / 10.0)
                        f_23 = (10 - f23) * f13 / 100.0
                        f_13 = f23 * (10 - f13) / 100.0
                        if t_intro_23 < t_intro_13:
                            f_23 += (10 - f23) * (10 - f13) / 100.0
                        else:
                            f_13 += (10 - f23) * (10 - f13) / 100.0
                        f_12 = 1 - f_23 - f_13
                        file_name = str_t12 + '_' + str_t123 + '_' + str_t_intro_23 + '_' + str_t_intro_13 + '_' + str(
                            '%.2f' % f_12) + '_' + str(f_13) + '_' + str(f_23) + '_' + '1323' + '_' + str(i)
                        OutputFile.write('seq_' + file_name + '\n')
                        TreeOutputFile.write('tree_' + file_name + '\n')
                        labels.write(str('%.2f' % f_12) + '	' + str(f_13) + '	' + str(f_23) + '\n')
                        if t_intro_23 <= t_intro_13:
                            os.system(
                                "%s %s -ej %s 2 1 -ej %s 3 1 -ej 3 4 1 -es %s 3 %s -ej %s 5 2 -es %s 3 %s -ej %s 6 1 %s -T | tail -n +4 | grep -v // > tree_%s" % (
                                    PATHTOMS, Indv_com, str_t12, str_t123, str_t_intro_23, str_f23, str_t_intro_23,
                                    str_t_intro_13, str_f13, str_t_intro_13, Rec_com, file_name))
                        else:
                            os.system(
                                "%s %s -ej %s 2 1 -ej %s 3 1 -ej 3 4 1 -es %s 3 %s -ej %s 5 1 -es %s 3 %s -ej %s 6 2 %s -T | tail -n +4 | grep -v // > tree_%s" % (
                                    PATHTOMS, Indv_com, str_t12, str_t123, str_t_intro_13, str_f13, str_t_intro_13,
                                    str_t_intro_23, str_f23, str_t_intro_23, Rec_com, file_name))

                        partitions = commands.getoutput('wc -l tree_' + file_name)
                        os.system("%s -q -mHKY %s -p %s < tree_%s > seq_%s" % (
                            PATHTOSEQGEN, Sub_com, partitions, file_name, file_name))

# 2019-11-04