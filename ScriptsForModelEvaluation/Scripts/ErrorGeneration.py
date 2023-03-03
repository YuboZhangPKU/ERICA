# _*_ coding: UTF-8 _*_

# Version information START --------------------------------------------------
VERSION_INFO = \
    """
    Author: ZHANG YUBO

    Version-01:
        2021-10-21 randomly introduce sequence error and missing data in geno file

    """
# Version information END ----------------------------------------------------

import argparse
import random
import gzip
###############################################################################
# read parameters
###############################################################################
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="introduce sequence error in geno")
    parser.add_argument("-i", "--Input",
                        help="Path to input geno file", required=True)
    parser.add_argument("-o", "--Output",
                        help="Output geno file", required=True)
    parser.add_argument("-w", "--WindowSize",
                        help="Length for each window", default="5000")
    parser.add_argument("-r", "--Ratio",
                        help="Error rate", default="0.01")
    parser.add_argument("-t", "--Type",
                        help="Error Type", choices={"Error", "Gap"}, default="Error")    

    # -------------------------------------------------------------------->>>>>
    # load the paramters
    # -------------------------------------------------------------------->>>>>
    ARGS = parser.parse_args()

    if ARGS.Input.endswith(".gz"):
        InputFile = gzip.open(ARGS.Input, "rt")
    else:
        InputFile = open(ARGS.Input, "r")

    OutputFile = ARGS.Output
    WindowSize = int(ARGS.WindowSize)
    Ratio = float(ARGS.Ratio)
    Count = int(Ratio * WindowSize)
    Type = ARGS.Type

    
    IndexList = list(range(WindowSize))

    lines = InputFile.readlines()
    with open(OutputFile, 'w') as FileOut:
        FileOut.write(lines[0])
        N = len(lines[0].strip().split())-2
        if Type =='Error':
            for i in range(1, len(lines), WindowSize):
                TempList = lines[i: i+WindowSize]
                for j in range(len(TempList)):
                    TempList[j] = TempList[j].strip().split()

                for j in range(N):
                    random.shuffle(IndexList)
                    for k in IndexList[0:Count]:
                        #print(k, end = ' ')
                        BaseList = ['A', 'T', 'C', 'G']
                        BaseList.remove(TempList[k][j+2])
                        TempList[k][j+2] = random.choice(BaseList)
                    #print('')
                for j in range(len(TempList)):
                    TempList[j] = '\t'.join(TempList[j])
                    FileOut.write(TempList[j]+'\n')    

        elif Type =='Gap':
            for i in range(1, len(lines), WindowSize):
                TempList = lines[i: i+WindowSize]
                for j in range(len(TempList)):
                    TempList[j] = TempList[j].strip().split()

                for j in range(N):
                    random.shuffle(IndexList)
                    for k in IndexList[0:Count]:
                        TempList[k][j+2] = 'N'

                for j in range(len(TempList)):
                    TempList[j] = '\t'.join(TempList[j])
                    FileOut.write(TempList[j]+'\n')  


# 2021-10-21