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


###############################################################################
# read parameters
###############################################################################
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="geno2vcf")
    parser.add_argument("-i", "--Input",
                        help="Input file name")
    parser.add_argument("-o", "--Output",
                        help="Output file name")
    parser.add_argument("-f", "--format",
                        help="Format to output",
                        choices=("haplo", "diplo", "homo"), required=True)


    # -------------------------------------------------------------------->>>>>
    # load the paramters
    # -------------------------------------------------------------------->>>>>
    ARGS = parser.parse_args()
    Input = ARGS.Input
    Output = ARGS.Output
    Format = ARGS.format

    if Input.endswith(".gz"):
        InputFile = gzip.open(Input, "rt")
    else:
        InputFile = open(Input, "r")

    with open(Output, 'w') as fout:
        if Format == "haplo":
            head = InputFile.readline().strip().split()
            HEAD = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + '\t'.join(head[2:]) + '\n'
            fout.write(HEAD)

            while True:
                line = InputFile.readline()
                if not line:
                    break
                line = line.strip().split()
                Ref = line[2]
                AllBase = list(set(line[2:]))
                if len(AllBase) > 1:
                    #print(AllBase)
                    AllBase.remove(Ref)
                    AllBase.insert(0, Ref)
                    Alt = ','.join(AllBase[1:])
                else:
                    Alt = '.'
                fout.write("%s\t%s\t.\t%s\t%s\t.\tPASS\t.\tGT" % (line[0], line[1], Ref, Alt))
                if Alt == '.':
                    fout.write("\t0"*(len(line)-2) + '\n')
                else:
                    for i in range(2, len(line)):
                        fout.write("\t" + str(AllBase.index(line[i])))
                    fout.write("\n")

        elif Format == "diplo":
            head = InputFile.readline().strip().split()
            newhead = []
            for i in range(2, len(head), 2):
                newhead.append(head[i])
            HEAD = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + '\t'.join(newhead) + '\n'
            fout.write(HEAD)

            while True:
                line = InputFile.readline()
                if not line:
                    break
                line = line.strip().split()
                Ref = line[2]
                AllBase = list(set(line[2:]))
                if len(AllBase) > 1:
                    #print(AllBase)
                    AllBase.remove(Ref)
                    AllBase.insert(0, Ref)
                    Alt = ','.join(AllBase[1:])
                else:
                    Alt = '.'
                fout.write("%s\t%s\t.\t%s\t%s\t.\tPASS\t.\tGT" % (line[0], line[1], Ref, Alt))
                if Alt == '.':
                    fout.write("\t0|0"* int((len(line)-2)/2) + '\n')
                else:
                    for i in range(2, len(line),2 ):
                        fout.write("\t" + str(AllBase.index(line[i])) + '|' +str(AllBase.index(line[i+1])))
                    fout.write("\n")                   
    
        elif Format == "homo":
            head = InputFile.readline().strip().split()
            HEAD = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + '\t'.join(head[2:]) + '\n'
            fout.write(HEAD)

            while True:
                line = InputFile.readline()
                if not line:
                    break
                line = line.strip().split()
                Ref = line[2]
                AllBase = list(set(line[2:]))
                if len(AllBase) > 1:
                    #print(AllBase)
                    AllBase.remove(Ref)
                    AllBase.insert(0, Ref)
                    Alt = ','.join(AllBase[1:])
                else:
                    Alt = '.'
                fout.write("%s\t%s\t.\t%s\t%s\t.\tPASS\t.\tGT" % (line[0], line[1], Ref, Alt))
                if Alt == '.':
                    fout.write("\t0|0"*(len(line)-2) + '\n')
                else:
                    for i in range(2, len(line)):
                        fout.write("\t" + str(AllBase.index(line[i])) + '|' + str(AllBase.index(line[i])))
                    fout.write("\n")