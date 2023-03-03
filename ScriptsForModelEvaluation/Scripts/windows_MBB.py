#!/usr/bin/env python
import sys,getopt,random,commands
try:
	opts, args = getopt.getopt(sys.argv[1:],"hi:o:w:b:",["ifile=","ofile="])
except getopt.GetoptError:
	print 'windows_MBB.py -i <input csv file> -o <output> -w <window size> -b  <block size>'
	sys.exit(2)
for opt, arg in opts:
	if opt == '-h':
		print 'windows_MBB.py -i <input csv file> -o <output> -w <window size> -b  <block size>'
		sys.exit()
	elif opt in ("-i", "--ifile"):
		input = arg
	elif opt in ("-o", "--ofile"):
		output = arg
	elif opt in ("-w"):
		window_size = arg
	elif opt in ("-b"):
		block_size = arg

fileout_1 = open(output+'_all_result.csv','w')
fileout_2 = open(output+'_positive_result.csv','w')
fileout_1.write('scaffold,start,D_mean,D_sd,D_pvalue,fd_mean,fd_sd,fd_pvalue\n')
fileout_2.write('scaffold,start,D_mean,D_sd,D_pvalue,fd_mean,fd_sd,fd_pvalue\n')
filein = open(input,'r')

rows = int(window_size)/int(block_size)
lines = filein.readlines()
for i in range(1,len(lines)):
	lines[i] = lines[i].split(',')
i = 1
stat = int(lines[i][1])
end =  stat + int(window_size)
flag = False

while True:
	#print(stat)
	#print(lines[i][1])
	#print('\n')
	temp_file = open('temp.csv','w') 
	for j in range(i,i+rows):
		if j >= len(lines):
			break
		if int(lines[j][1]) >= end:
			j = j-1
			break
		if lines[j][0] != lines[i][0]:
			flag = True
			j = j-1
			break
		else:
			if lines[j][8] == 'nan' or lines[j][8] == '-inf':
				D = '0'
			else:
				D = lines[j][8]
			if lines[j][9] == 'nan' or lines[j][9] == '-inf':
				fd = '0'
			else:
				fd = lines[j][9]

			temp_file.write(lines[j][1]+','+D+','+fd+'\n')
		
	temp_file.close()
	
	result = commands.getoutput("Rscript moving_block_bootstrap.R")
	#print(result)
	result = result.split()
	
	write_result = lines[i][0]+','+str(stat)
	for value in result:
		write_result += ','
		write_result += value
	write_result += '\n'
	fileout_1.write(write_result)
	if float(result[2]) < 0.05 and float(result[5]) < 0.05:
		fileout_2.write(write_result)

	i = j+1
	if flag:
		stat = int(lines[i][1])
		end =  stat + int(window_size)
		flag = False
	else:
		stat += int(window_size)
		end += int(window_size)
	if i >= len(lines):
		break

filein.close()	
fileout_1.close()
fileout_2.close()