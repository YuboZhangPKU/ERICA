# _*_ coding: UTF-8 _*_
# Version information START --------------------------------------------------
VERSION_INFO = \
    """
    Author: ZHANG YUBO

    Version-01:
        2020-01  Post-processing and visualization of the ERICA results

    Version-02:
        2020-08  Processing both four-taxon and five-taxon results         
    """
# Version information END ----------------------------------------------------


import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas
from plotnine import *


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Summary and visualization of ERICA results")
    parser.add_argument("-i", "--InputFile",
                        help="Input file name of ERICA results", required=True)
    parser.add_argument("-o", "--OutputFile",
                        help="Output prefix", required=True)
    parser.add_argument("-p", "--Population",
                        help="Number of populations", choices=("4", "5"), required=True)
    parser.add_argument("-w", "--WindowSize",
                        help="Window size for outputting and plotting, default=50000", default="50000")
    parser.add_argument("-r", "--Region",
                        help="1-based indexes for analyzed regions, default=1:-1 (the full region)", default='1:-1')
    parser.add_argument("-d", "--DistanceToZero",
                        help="The threshold used for classification, probabilities greater than the threshold will be recorded", default='0.40')
    parser.add_argument("-m", "--MaxValue",
                        help="Plot the highest supporting topologies along chromosome", default='True', choices=("True", "False"))
    parser.add_argument("-l", "--Line",
                        help="Line plot of each topology along chromosome", default='True', choices=("True", "False"))
    parser.add_argument("-a", "--Area",
                        help="Area plot of each topology along chromosome", default='True', choices=("True", "False"))
    parser.add_argument("-c", "--Chr",
                        help="Name for the chromosome, optional", required=False)


    # -------------------------------------------------------------------->>>>>
    # load the paramters
    # -------------------------------------------------------------------->>>>>
    ARGS = parser.parse_args()
    InputFile = ARGS.InputFile
    OutputFile = ARGS.OutputFile
    PopulationCount = ARGS.Population
    raw_window_size = 5000
    width = int(ARGS.WindowSize) // raw_window_size


    region=ARGS.Region.split(':')
    with open (ARGS.InputFile,'r') as filein:
        raw_data = filein.readlines()
    for i in [0,1]:
        if '-' in region[i]:
            region[i] = int(region[i]) + len(raw_data) + 1
        else:
            region[i] = int(region[i])
    raw_data = raw_data[region[0] - 1 : region[1]]
    index_list = range(((region[0] - 1) * raw_window_size + 1), (((region[1] - 1) // width + 1) * raw_window_size * width + 1), raw_window_size * width)


    DistanceToZero = float(ARGS.DistanceToZero)
    MaxValueFlag = ARGS.MaxValue
    LineFlag=ARGS.Line
    AreaFlag=ARGS.Area
    if ARGS.Chr:
        chr_name = ARGS.Chr

    
    #calculating mean value for each window
    mean_data = []
    if PopulationCount == "4":
        n = 3
        TopoName = ["A", "B", "C"]
        color_list = ["#ECA257","#B3B2B2","#1766A0"]
    elif PopulationCount == "5":
        n = 15
        TopoName = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O"]
        color_list = ["#06b8b9", "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928", "#b7b7b3", "#1e2022"]


    raw_data_count = 0
    window_sum = [0] * n
    for raw_data_line in raw_data:
            raw_data_line = raw_data_line.split()
            raw_data_count += 1
            if raw_data_count < width:
                for index in range(n):
                    window_sum[index] += float(raw_data_line[index])
            else:
                temp_list = []
                for index in range(n):
                    window_sum[index] += float(raw_data_line[index])
                    temp_list.append(window_sum[index]/width)
                mean_data.append(temp_list)
                raw_data_count = 0
                window_sum = [0] * n
    if raw_data_count != 0:
        temp_list = []
        for index in range(n):
            temp_list.append(window_sum[index] / raw_data_count)
        mean_data.append(temp_list)


    mean_data_df = pandas.DataFrame(mean_data)
    mean_data_df.columns = TopoName
    mean_data_df["Index"] = pandas.Series(index_list)
    if ARGS.Chr:
        mean_data_df["Chr"] = pandas.Series([chr_name]*len(index_list))
        order = ['Chr', 'Index'] + TopoName
        mean_data_df = mean_data_df[order]
    else:
        order = ['Index'] + TopoName
        mean_data_df = mean_data_df[order]


    # topo C minus topo B, compared with ABBA-BABA test
    if PopulationCount == "4":
        mean_data_df["CminusB"] = mean_data_df["C"] - mean_data_df["B"]


    max_value = []
    max_class = []
    max_data = []
    for mean_value_list in mean_data:
        max_value.append(max(mean_value_list))
        max_class.append(chr(65+mean_value_list.index(max(mean_value_list))))
        temp_list = [0] * n
        temp_list[mean_value_list.index(max(mean_value_list))] = max(mean_value_list)
        max_data.append(temp_list)
    mean_data_df["Max Value"] = pandas.Series(max_value)
    mean_data_df["Max Class"] = pandas.Series(max_class)
    distance_class = []
    for i in mean_data_df.index:
        temp_class = ''
        for j in TopoName:
            if mean_data_df[j][i] > DistanceToZero:
                temp_class += j
        if temp_class == '':
            temp_class = '?'
        distance_class.append(temp_class)
    mean_data_df["Distance Class"] = pandas.Series(distance_class)
    mean_data_df.to_csv(OutputFile+'.csv', index = False)


    if MaxValueFlag == 'True':
        max_data_df = pandas.DataFrame(max_data)
        p = ggplot(data = max_data_df)
        for i in range(n):
            p += geom_area(aes(x = max_data_df.index, y = max_data_df[i]), color = color_list[i],fill = color_list[i])
        p += scale_y_continuous(limits= [0,1])
        p += xlab("")
        p += ylab("")
        p += theme_classic()
        pdf_height = 15
        pdf_width = 40
        plot_file_name = OutputFile + '_MaxValue.pdf'
        p.save(plot_file_name, height=pdf_height, width=pdf_width,limitsize=False )


    if AreaFlag ==  'True':
        p = ggplot(data = mean_data_df)
        for i in range(n - 1, -1, -1):
            temp_list = mean_data_df["A"].tolist()
            temp_serie = pandas.Series(temp_list)
            for j in range(1,i+1):
                temp_serie += mean_data_df[TopoName[j]]
            p += geom_area(aes(x = mean_data_df.index, y = temp_serie), color = color_list[i],fill = color_list[i])
        p += scale_y_continuous(limits= [0,1])
        p += xlab("")
        p += ylab("")
        p += theme_classic()
        pdf_height = 15
        pdf_width = 40
        plot_file_name = OutputFile + '_Area.pdf'
        p.save(plot_file_name, height=pdf_height, width=pdf_width,limitsize=False )


    if LineFlag ==  'True':
        p = ggplot(data = mean_data_df)
        for i in range(n):
            p += geom_line(aes(x = mean_data_df.index, y = mean_data_df[TopoName[i]]), color = color_list[i], size = 2)
        p += scale_y_continuous(limits= [0,1])
        p += xlab("")
        p += ylab("")
        p += theme_classic()
        pdf_height = 15
        pdf_width = 40
        plot_file_name = OutputFile + '_Line.pdf'
        p.save(plot_file_name, height=pdf_height, width=pdf_width,limitsize=False )

