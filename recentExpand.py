import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("-i", "--workdir", required=True)
parser.add_argument("-f", "--input_file", required=False, default='/pattern.fa')
parser.add_argument("-n", "--cluster_number", required=True)
args = parser.parse_args()

workdir = args.workdir
file = workdir + args.input_file
out_file = workdir + '/match_in_window.xls'

sample_seq = {}

with open(file, 'r') as f:
    while True:
        line = f.readline()[:-1]
        if not line:
            break
        name = line[1:].split('_')
        seq = f.readline()[:-1]
        print(name)
        if name[0].startswith('CHM'):
            if name[0] not in sample_seq.keys():
                start = int(name[3].split('-')[0])
                sample_seq[name[0]] = [[int(name[3].split('-')[0]),int(name[3].split('-')[1]),seq]]
            else:
                sample_seq[name[0]].append([int(name[3].split('-')[0]),int(name[3].split('-')[1]),seq])
        else:
            if name[0]+name[1] not in sample_seq.keys():
                start = int(name[3].split('-')[0])
                sample_seq[name[0]+name[1]] = [[int(name[3].split('-')[0]),int(name[3].split('-')[1]),seq]]
            else:
                sample_seq[name[0]+name[1]].append([int(name[3].split('-')[0]),int(name[3].split('-')[1]),seq])

sort_sample_seq = {}
for i in sample_seq.keys():
    sort_sample_seq[i] = sorted(sample_seq[i],key=lambda x:x[0])

# 划分window,1unit滑动 10个一组检查其中完全一致的个数
split_sample_seq = {}
sample_window_of_match = {}
window_size = 10
for i in sort_sample_seq.keys():
    split_sample_seq[i] = []
    sample_window_of_match[i] = []

    for j in range(len(sort_sample_seq[i]) - window_size):

        split_sample_seq[i].append(sort_sample_seq[i][j:j+window_size])

    for j in split_sample_seq[i]:
        start = j[0][0]
        end = j[0][1]
        match_number = 0
        for k in j:
            for l in j:
                if k[0] == l[0]:
                    continue
                if k[2] == l[2]:
                    match_number += 1
        sample_window_of_match[i].append([start,end,match_number])


out_file = open(out_file,'w')
for i in sample_window_of_match.keys():
    sample = i
    for j in sample_window_of_match[i]:
        start = j[0]
        end = j[1]
        number = j[2]
        out_file.write(sample + '\t' + str(start)+'\t'+str(end)+'\t'+str(number)+'\n')
out_file.close()

