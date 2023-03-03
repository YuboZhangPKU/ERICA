#!/usr/bin/env python
import argparse
import os

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calculate Dfoil using sliding windows")
    parser.add_argument("-i", "--InputFile",
                        help="Input topo file", required=True)
    parser.add_argument("-o", "--OutputFile",
                        help="Output ", required=True)
    parser.add_argument("-w", "--WindowSize",
                        help="Window size, default=5000", default="5000")
    parser.add_argument("-P1", "--pop1", help="Sample name", required=True)
    parser.add_argument("-P2", "--pop2", help="Sample name", required=True)
    parser.add_argument("-P3", "--pop3", help="Sample name", required=True)
    parser.add_argument("-P4", "--pop4", help="Sample name", required=True)
    parser.add_argument("-O", "--outgroup", help="Sample name", required=True)
    # -------------------------------------------------------------------->>>>>
    # load the paramters
    # -------------------------------------------------------------------->>>>>
    ARGS = parser.parse_args()

    InputFile=ARGS.InputFile
    OutputFile=ARGS.OutputFile
    WindowSize = int(ARGS.WindowSize)
    pop1=ARGS.pop1
    pop2=ARGS.pop2
    pop3=ARGS.pop3
    pop4=ARGS.pop4
    outgroup=ARGS.outgroup

    start = 1

    pop_list = [pop1,pop2,pop3,pop4,outgroup]
    seq_dict = {}

    with open(InputFile,'r') as file_in:
        while True:
            temp_name = file_in.readline()
            if not temp_name:
                break
            temp_name = temp_name.strip()[1:]
            temp_seq = file_in.readline()
            if temp_name in pop_list:
                seq_dict[temp_name]=temp_seq.strip()
    temp_seq =''

    sequence_len = len(seq_dict[pop1])

    file_out = open(OutputFile,'w')

    while True:
        if start+WindowSize-1 > sequence_len:
            break
        temp_file_out = open('temp.fa','w')
        for i in range(5):
            temp_file_out.write('>' + str(i+1) + '\n')
            temp_file_out.write(seq_dict[pop_list[i]][start-1:start-1+WindowSize])
            temp_file_out.write('\n')
        temp_file_out.close()
        os.system('python /home/zhangyubo/Simulation/dfoil-master/fasta2dfoil.py temp.fa -o temp.counts --names 1,2,3,4,5')
        os.system('python /home/zhangyubo/Simulation/dfoil-master/dfoil.py --infile temp.counts --out temp.dfoil --mode dfoil')
        temp_file_in = open('temp.dfoil', 'r')
        lines = temp_file_in.readlines()
        temp_file_in.close()
        try:
            file_out.write(lines[1])
        except IndexError:
            file_out.write('nan\n')

        for file in ['temp.counts', 'temp.fa', 'temp.dfoil']:
            os.system('rm -rf ' + file)

        start += WindowSize

