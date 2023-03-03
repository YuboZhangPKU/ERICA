# _*_ coding: UTF-8 _*_

# Version information START --------------------------------------------------
VERSION_INFO = \
    """
    Author: ZHANG YUBO

    Version-01:
        2019-08-18 simulate asymmetric five-taxon data set

    Version-02:
        2019-11-04 simulate asymmetric five-taxon data set, with variable mutation rates, recombination rates and number of haplotype

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
    parser = argparse.ArgumentParser(description="Simulate asymmetric five-taxon data set")
    parser.add_argument("-o", "--Output",
                        help="Path to output directory")
    parser.add_argument("--RepNum",
                        help="Repeat all of the 518 evolutionary scenarios 10 times for training dataset and 1 times for test dataset",
                        default="10")
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

    RepNum = int(ARGS.RepNum)
    RecRate = float(ARGS.RecRate)
    SubRate = float(ARGS.SubRate)
    NumPerTaxon = int(ARGS.NumPerTaxon)
    Length = int(ARGS.Length)
    PATHTOMS = ARGS.PathToMs
    PATHTOSEQGEN = ARGS.PathToSeqgen

    OutputFile = open('file_name', 'w')
    labels = open('labels', 'w')
    TreeOutputFile = open('tree_file_name', 'w')

    Indv_com = str(5 * NumPerTaxon) + ' 1 -I 5 ' + (str(NumPerTaxon) + ' ') * 5
    if RecRate != 0:
        Rec_com = '-r %s %s' % (int(RecRate * Length), Length)
    else:
        Rec_com = ''
    Sub_com = '-l %s -s %s' % (Length, SubRate)

    # part1 topo((((1,2),3),4),5)
    t12_list = [2, 5, 8, 11, 14, 17]
    for i in range(1, RepNum + 1):
        for t12 in t12_list:
            for t123 in range(t12 + 2, 23, 3):
                for t1234 in range(t123 + 2, 28, 3):
                    file_name = "%s_%s_%s_%s" % ((t12/10.0), (t123/10.0), (t1234/10.0), i)
                    OutputFile.write('seq_%s\n' % file_name)
                    TreeOutputFile.write('tree_%s\n' % file_name)
                    labels.write('1 0 0 0 0 0 0 0 0 0 0 0 0 0 0' + '\n')
                    os.system("%s %s -ej %s 2 1 -ej %s 3 1 -ej %s 4 1 -ej 3 5 1 %s -T | tail -n +4 | grep -v // > tree_%s" % (
                        PATHTOMS, Indv_com, (t12/10.0), (t123/10.0), (t1234/10.0), Rec_com, file_name))
                    partitions = commands.getoutput('wc -l tree_' + file_name)
                    os.system(
                        "%s -q -mHKY %s -p %s < tree_%s > seq_%s" % (
                        PATHTOSEQGEN, Sub_com, partitions, file_name, file_name))

    # part2 topo((((1,2),3),4),5) one introgression
    intro_list = ['intro13', 'intro31', 'intro23', 'intro32', 'intro14', 'intro41', 'intro24', 'intro42', 'intro124',
                  'intro412', 'intro34', 'intro43']
    com_list = ['31', '13', '32', '23', '41', '14', '42', '24', '41', '14', '43', '34']
    topo_list = ['2', '2', '3', '3', '11', '15', '12', '14', '10', '10', '13', '13']
    for j in range(0, 12):
        for i in range(1, RepNum + 1):
            t12 = random.randint(4, 15)
            t123 = random.randint(t12 + 2, 21)
            t1234 = random.randint(t123 + 2, 27)
            if j < 8:
                t_intro = random.randint(1, t12 - 1)
            elif j < 10:
                t_intro = random.randint(t12 + 1, t123 - 1)
            else:
                t_intro = random.randint(1, t123 - 1)

            for f in range(0, 10):
                file_name = "%s_%s_%s_%s_%s_%s_%s" % ((t12/10.0), (t123/10.0), (t1234/10.0), (t_intro/10.0), (f/10.0), intro_list[j], i)
                OutputFile.write('seq_%s\n' % file_name)
                TreeOutputFile.write('tree_%s\n' % file_name)
                Intro_com = '-es %s %s %s -ej %s 6 %s ' % ((t_intro/10.0), com_list[j][0], (f/10.0), (t_intro/10.0), com_list[j][1])
                os.system(
                    "%s %s -ej %s 2 1 -ej %s 3 1 -ej %s 4 1 -ej 3 5 1 %s %s -T | tail -n +4 | grep -v // > tree_%s" % (
                        PATHTOMS, Indv_com, (t12/10.0), (t123/10.0), (t1234/10.0), Intro_com , Rec_com, file_name))
                partitions = commands.getoutput('wc -l tree_' + file_name)
                os.system(
                    "%s -q -mHKY %s -p %s < tree_%s > seq_%s" % (
                        PATHTOSEQGEN, Sub_com, partitions, file_name, file_name))
                write_labels = '%s %s%s%s\n' % ((f/10.0), ('0 ' * (int(topo_list[j])-2)), (1-(f/10.0)), (' 0' * (15 - int(topo_list[j]))))
                labels.write(write_labels)

    # part3 topo((((1,2),3),4),5) two introgression
    intro_list = ['intro13', 'intro31', 'intro23', 'intro32', 'intro14', 'intro41', 'intro24', 'intro42', 'intro124',
                  'intro412', 'intro34', 'intro43']
    com_list = ['31', '13', '32', '23', '41', '14', '42', '24', '41', '14', '43', '34']
    topo_list = ['2', '2', '3', '3', '11', '15', '12', '14', '10', '10', '13', '13']
    two_intro_topo_list = [['2', '2', '2', '2', '7', '7', '14', '14', '2', '2', '2', '2'],
                           ['2', '2', '2', '2', '14', '2', '14', '14', '14', '14', '7', '7'],
                           ['3', '3', '3', '3', '15', '15', '4', '4', '3', '3', '3', '3'],
                           ['3', '3', '3', '3', '15', '15', '15', '3', '15', '15', '4', '4'],
                           ['8', '8', '15', '15', '11', '15', '11', '8', '11', '11', '11', '11'],
                           ['15', '15', '15', '15', '11', '15', '11', '11', '11', '11', '8', '8'],
                           ['14', '14', '5', '5', '12', '5', '12', '14', '12', '12', '12', '12'],
                           ['14', '14', '14', '14', '12', '12', '12', '14', '12', '12', '5', '5'],
                           ['2', '14', '3', '15', '11', '11', '12', '12', '10', '10', '13', '13'],
                           ['2', '14', '3', '15', '11', '11', '12', '12', '10', '10', '13', '13'],
                           ['9', '9', '6', '6', '13', '6', '13', '9', '13', '13', '13', '13'],
                           ['13', '13', '13', '13', '9', '9', '6', '6', '13', '13', '13', '13']]

    for j in range(0, 12):
        for k in range(0, 12):
            for i in range(1, 2 * RepNum + 1):
                t12 = random.randint(10, 15)
                t123 = random.randint(t12 + 3, 22)
                t1234 = random.randint(t123 + 2, 27)
                if j < 8:
                    t_intro_1 = random.randint(5, t12 - 1)
                elif j < 10:
                    t_intro_1 = random.randint(t12 + 1, t123 - 2)
                else:
                    t_intro_1 = random.randint(5, t123 - 2)
                if k < 8:
                    t_intro_2 = random.randint(1, min(t_intro_1, t12) - 1)
                elif k < 10:
                    t_intro_2 = random.randint(max(t12 + 1, t_intro_1 + 1), t123 - 1)
                else:
                    t_intro_2 = random.randint(1, min(t_intro_1, t123) - 1)
                f_1 = random.randint(0, 10)
                f_2 = random.randint(0, 10)

                file_name = "%s_%s_%s_%s_%s_%s_%s_%s_%s_%s" % ((t12/10.0), (t123/10.0), (t1234/10.0), (t_intro_1/10.0), (f_1/10.0), intro_list[j], (t_intro_2/10.0), (f_2/10.0), intro_list[k],i)
                OutputFile.write('seq_%s\n' % file_name)
                TreeOutputFile.write('tree_%s\n' % file_name)
                if t_intro_1 < t_intro_2:
                    Intro_com = '-es %s %s %s -ej %s 6 %s ' % (
                    (t_intro_1 / 10.0), com_list[j][0], (f_1 / 10.0), (t_intro_1 / 10.0), com_list[j][1])
                    Intro_com += '-es %s %s %s -ej %s 7 %s ' % (
                    (t_intro_2 / 10.0), com_list[k][0], (f_2 / 10.0), (t_intro_2 / 10.0), com_list[k][1])
                else:
                    Intro_com = '-es %s %s %s -ej %s 6 %s ' % (
                    (t_intro_2 / 10.0), com_list[k][0], (f_2 / 10.0), (t_intro_2 / 10.0), com_list[k][1])
                    Intro_com += '-es %s %s %s -ej %s 7 %s ' % (
                    (t_intro_1 / 10.0), com_list[j][0], (f_1 / 10.0), (t_intro_1 / 10.0), com_list[j][1])
                os.system(
                    "%s %s -ej %s 2 1 -ej %s 3 1 -ej %s 4 1 -ej 3 5 1 %s %s -T | tail -n +4 | grep -v // > tree_%s" % (
                        PATHTOMS, Indv_com, (t12 / 10.0), (t123 / 10.0), (t1234 / 10.0), Intro_com, Rec_com,
                        file_name))
                partitions = commands.getoutput('wc -l tree_' + file_name)
                os.system(
                    "%s -q -mHKY %s -p %s < tree_%s > seq_%s" % (
                        PATHTOSEQGEN, Sub_com, partitions, file_name, file_name))

                write_labels = [0] * 16
                write_labels[1] += (f_1 / 10.0) * (f_2 / 10.0)
                write_labels[int(topo_list[k])] += (f_1 / 10.0) * (1 - (f_2 / 10.0))
                write_labels[int(topo_list[j])] += (1 - (f_1 / 10.0)) * (f_2 / 10.0)
                write_labels[int(two_intro_topo_list[k][j])] += (1 - (f_1 / 10.0)) * (1 - (f_2 / 10.0))
                write_label = ''
                for writelabel in write_labels[1:]:
                    write_label += str(writelabel) + ' '
                labels.write(write_label[:-1] + '\n')
# 2019-12-01