#!/usr/bin/env python 
import sys
import json

sample_name = sys.argv[1]
file = sys.argv[2]

with open(file, 'r') as data_file:
    data = json.load(data_file)

with open(sample_name+"_duplication-rate.txt", 'w') as dup_file:
    dup_file.write("sample_name\tn_reads\tdup_reads\n")
    for point in data[0]['Fragment']['duplicate_saturation']:
        line = "\t".join([
            sample_name,
            str(point[0]),
            str(point[1])
        ])
        dup_file.write(line+"\n")
        