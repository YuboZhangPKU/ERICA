# _*_ coding: UTF-8 _*_
import argparse
import os
import commands
from multiprocessing import Process
import os
import time
import re
import subprocess

###############################################################################
# global variables
PATH_GPHOCS = '/usr/share/G-PhoCS/bin/G-PhoCS'

###############################################################################
def run_gphocs(input_ctl, outputname, mig_band):
    ctl_file_name =  '%s_mig_%s2%s.ctl' % (outputname,mig_band['source'],mig_band['target'])
    with open(input_ctl,'r') as fin:
        lines = fin.readlines()
    with open(ctl_file_name, 'w') as fout:
        for line in lines:
            try:
                note = line.split()[0]
                if note == 'trace-file':
                    fout.write('\ttrace-file\t\t%s_mig_%s2%s.log\n'% (outputname,mig_band['source'],mig_band['target']) )
                elif note == 'MIG-BANDS-START':
                    fout.write(line)
                    fout.write('\tBAND-START\n\tsource\t%s\n\ttarget\t%s\n\tBAND-END\n' %(mig_band['source'],mig_band['target']))
                else:
                    fout.write(line)
            except:
                fout.write(line)

    popen = subprocess.Popen('%s %s' % (PATH_GPHOCS, ctl_file_name),
                             stdin=subprocess.PIPE,
                             shell=True)
    popen.communicate("\n")

###############################################################################
# read parameters
###############################################################################
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="run migration model of gphocs")
    parser.add_argument("-i", "--Input",
                        help="Input control file without migration bands")
    parser.add_argument("-o", "--Output",
                        help="Output file name")
    parser.add_argument("-p", "--Population",
                        help="Population name, split with comma")
    parser.add_argument("-T", "--Thread", default="4",
                        help="Thread for multiprocessing")
    # -------------------------------------------------------------------->>>>>
    # load the paramters
    # -------------------------------------------------------------------->>>>>
    ARGS = parser.parse_args()
    full_cmd = """python run_gphocs.py\\
    -i {Input} \\
    -o {Output} \\
    -p {Population} \\
    -T {Thread}
    """.format(
        Input=ARGS.Input,
        Output=ARGS.Output,
        Population=ARGS.Population,
        Thread=ARGS.Thread
    )

    Input = ARGS.Input
    Output=ARGS.Output
    Population=ARGS.Population.split(',')
    Thread = int(ARGS.Thread)
    migration_model_list =[]
    for i in Population:
        for j in Population:
            if i != j:
                migration_model_list.append({'source':i,'target':j})
    print(migration_model_list)
    index = 0
    flag = True
    while flag:
        processes = []
        for i in range(Thread):
            try:
                p = Process(target=run_gphocs, args=(Input, Output, migration_model_list[index]))
                print('Process will start.')
                p.start()
                processes.append(p)
                index += 1
            except:
                flag = False
                continue

        for p in processes:
            p.join()
        print('Process end.')


# 2020-1-3