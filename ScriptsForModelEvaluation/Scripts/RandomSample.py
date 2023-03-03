# _*_ coding: UTF-8 _*_

# Version information START --------------------------------------------------
VERSION_INFO = \
    """
    Author: ZHANG YUBO

    Version-01:
        2021-10-9 geno to vcf
    """
# Version information END ----------------------------------------------------

import argparse
import random
import gzip

###############################################################################
# read parameters
###############################################################################
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="geno2vcf")
    parser.add_argument("-i", "--Input",
                        help="Input file name")
    parser.add_argument("-o", "--Output",
                        help="Output file name")
    parser.add_argument("--Indv",
                        help="")
    parser.add_argument("-w", "--WindowSize",
                        help="Window size", default="5000")


    # -------------------------------------------------------------------->>>>>
    # load the paramters
    # -------------------------------------------------------------------->>>>>
    ARGS = parser.parse_args()
    Input = ARGS.Input
    Output = ARGS.Output
    Indv = ARGS.Indv.split(',')
    WindowSize = int(ARGS.WindowSize)

    if Input.endswith(".gz"):
        InputFile = gzip.open(Input, "rt")
    else:
        InputFile = open(Input, "r")

    with open(Output, 'w') as fout:
        HEAD = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tArchaic\n"
        fout.write(HEAD)

        head = InputFile.readline().split()
        IndvIndex = []
        for i in Indv:
            IndvIndex.append(head.index(i))

        SiteIndex = 0

        while True:
            line = InputFile.readline()
            if not line:
                break    
            if SiteIndex % WindowSize == 0:
                SampleID = random.choice(IndvIndex)
            #print(SampleID)
            line = line.split()
            fout.write("\t".join(line[:9]) + "\t")
            fout.write(line[SampleID]+ "\n")
            SiteIndex += 1
