#!/usr/bin/python3

# 2022-06-20 version 2

import sys
import os
import re
import csv
from bs4 import BeautifulSoup
import collections

# command verification 
args = sys.argv

if len(args) !=2:
    sys.exit("\nError\nUsage: python3 [program].py [input file name]\n")

# initial parameters
file_input = args[1]
file_prefix = os.path.splitext(os.path.basename(file_input))[0]
path_input = os.getcwd()

# file open
with open(file_input, "r") as f:
    if not 'xml' in f.readline():
        sys.exit("\nError\nInput file is not XML file.\nRequirement: blast output file with run command option '-outfmt 16' (Single-file BLAST XML2)\n")
    else:
        print("\nFiletype is xml. OK\n")
        blastresult = f.read()

# soup
print("\nReading file\n")

soup_1 = BeautifulSoup(blastresult, "lxml-xml")
result_list = soup_1.find_all('BlastOutput2')

# extract and file output
output_file1 = file_prefix + '_blasthittax.tsv'
output_file2 = file_prefix + '_blasthitsummary.tsv'

print("\nExtracting\n")

for result in result_list:
    query_acc = result.find('query-title').text
    query_len = result.find('query-len').text
    query_list = str.splitlines(query_acc) + str.splitlines(query_len)
    with open(output_file1, "a") as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(query_list)
    hit_acc = []
    hit_taxid = []
    hit_sciname = []
    hit_len = []
    for accession in result.find_all('accession'):
        hit_acc.append(accession.text)
    for taxid in result.find_all('taxid'):
        hit_taxid.append(taxid.text)
    for sciname in result.find_all('sciname'):
        hit_sciname.append(sciname.text)
    for length in result.find_all('len'):
        hit_len.append(length.text)
    cont_num = len(list(hit_acc))
    query_summary = str.splitlines(query_acc) + str.splitlines(str(cont_num))
    total = collections.Counter(hit_sciname)
    with open(output_file2, "a") as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(query_summary)
        for k, v in total.items():
            hit_summary = str.splitlines(str(k)) + str.splitlines(str(v))
            writer.writerow(hit_summary)
    for i in range(cont_num):
        hit_list = hit_acc[i:i+1] + hit_taxid[i:i+1] + hit_sciname[i:i+1] + hit_len[i:i+1]
        with open(output_file1, "a") as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(hit_list)

print("\nDone\n")
