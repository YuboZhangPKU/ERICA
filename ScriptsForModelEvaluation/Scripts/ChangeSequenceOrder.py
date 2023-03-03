# Version information START --------------------------------------------------
VERSION_INFO = \
    """
    Author: ZHANG YUBO

    Version-01:
        2019-08-18 change sequence order to generate alternative topology

    """
# Version information END ----------------------------------------------------

import argparse
import os
###############################################################################
# read parameters
###############################################################################
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Change sequence order")
    parser.add_argument("-i", "--Input",
                        help="Path to input directory")
    parser.add_argument("-o", "--Output",
                        help="Path to output directory")
    parser.add_argument("--Order",
                        help="New population order")
    parser.add_argument("--NumOfTaxon",
                        help="Number of taxon", default="4")
    parser.add_argument("--NumPerTaxon",
                        help="Number of samples per taxon", default="8")
    parser.add_argument("--NameFile",
                        help="Path to the file recording file names to be analysed", default="file_name")

    # -------------------------------------------------------------------->>>>>
    # load the paramters
    # -------------------------------------------------------------------->>>>>
    ARGS = parser.parse_args()
    InputDir = ARGS.Input
    OutputDir = ARGS.Output
    if not os.path.exists(OutputDir):
        os.mkdir(OutputDir)
    
    Order = ARGS.Order
    NumOfTaxon = int(ARGS.NumOfTaxon)
    NumPerTaxon = int(ARGS.NumPerTaxon)
    NameFile = ARGS.NameFile
    os.system("cp %s/%s %s/%s" % (InputDir, NameFile, OutputDir, NameFile))
    
    with open(InputDir+'/'+NameFile,'r') as NameFileInput:
        FileNames = NameFileInput.readlines()

    for FileName in FileNames:
        FileName = FileName.strip()
        with open(InputDir+'/'+FileName,'r') as TempFileInput:
            lines = TempFileInput.readlines()

        with open(OutputDir+'/'+FileName,'w') as TempFileOutput:
            TempFileOutput.write(lines[0])
            lines = lines[1:]
            for i in range(len(lines)):
                lines[i] = lines[i].split()
            lines.sort(key=lambda x: int(x[0]))

            if NumOfTaxon == 4:
                OriginOrder = '1234'
                for i in range(4):
                    for j in range(NumPerTaxon):
                        lines[j + (int(OriginOrder[i]) - 1) * NumPerTaxon][0] = str(j + (int(Order[i]) - 1) * NumPerTaxon + 1)
                for line in lines:
                    TempFileOutput.write(line[0] + '\t' + line[1] + '\n')

            elif NumOfTaxon == 5:
                OriginOrder = '12345'
                for i in range(5):
                    for j in range(NumPerTaxon):
                        lines[j + (int(OriginOrder[i]) - 1) * NumPerTaxon][0] = str(j + (int(Order[i]) - 1) * NumPerTaxon + 1)
                for line in lines:
                    TempFileOutput.write(line[0] + '\t' + line[1] + '\n')

# 2019-08-18
