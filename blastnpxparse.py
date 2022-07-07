#!/usr/bin/python3

# 2022-07-07 version 6.0: blastnpxparse.py

import sys
import os
#import re
import subprocess
import csv
from bs4 import BeautifulSoup
import collections

# command verification 
args = sys.argv

if len(args) !=2:
    sys.exit("\nError\nUsage: python3 [program].py [input file name]\n")

# initial parameters
file_input = args[1]
file_prefix_init = os.path.splitext(os.path.basename(file_input))[0]
file_conversion = '0'
path_input = os.getcwd()

# file open
with open(file_input, "r", encoding="utf-8") as f:
    headline = f.readline()
    if 'xml' in headline:
        print("\nInput file type is XML.\n\n\nReading file ...\n")
        file_prefix = file_prefix_init
        blastxml = f.read()
        if 'blastn' in blastxml:
            mode = 'blastn'
            print("\nblastn output file confirmed.\n")
        elif 'blastp' in blastxml:
            mode = 'blastp'
            print("\nblastp output file confirmed.\n")
        elif 'blastx' in blastxml:
            mode = 'blastx'
            print("\nblastx output file confirmed.\n")
    elif 'Blast4-archive' in headline:
        file_conversion = '1'
        print("\nInput file type is ASN.1.\n\n\nStarting file conversion to XML (Single-file BLAST XML2)\n")
        pass
    else:
        sys.exit("\nError\nInput is neither XML nor ASN.1 file.\nRequirement: blast[n, p, x] output file with run command option '-outfmt 11' (BLAST archive (ASN.1)) or '-outfmt 16' (Single-file BLAST XML2)\n")

if file_conversion == '1':
    file_conv = file_prefix_init + '_16' + os.path.splitext(os.path.basename(file_input))[1]
    file_prefix = file_prefix_init + '_16'
    cmd_list = ['blast_formatter','-archive',file_input,'-out',file_conv,'-outfmt',16]
    cmd = map(str, cmd_list)
    blastformatter_run = subprocess.run(cmd)
    print("\nblast_formatter finished.\n\n\nReading file ...\n")
    with open(file_conv, "r", encoding="utf-8") as f:
        blastxml = f.read()
        if 'blastn' in blastxml:
            mode = 'blastn'
            print("\nblastn output file confirmed.\n")
        elif 'blastp' in blastxml:
            mode = 'blastp'
            print("\nblastp output file confirmed.\n")
        elif 'blastx' in blastxml:
            mode = 'blastx'
            print("\nblastx output file confirmed.\n")
        else:
            sys.exit("\nError\nInput is not blast[n, p, x] file.\nRequirement: blast[n, p, x] output file with run command option '-outfmt 11' (BLAST archive (ASN.1)) or '-outfmt 16' (Single-file BLAST XML2)\n")
elif file_conversion == '0':
    pass

# soup
soup_1 = BeautifulSoup(blastxml, "lxml-xml")
result_queries = soup_1.find_all('BlastOutput2')

# extract and file output
output_file1 = file_prefix + '_' + mode + 'table.tsv'
output_file2 = file_prefix + '_' + mode + 'summary.tsv'

print("\nExtracting ...\n")

for result_query in result_queries:
    if result_query.find('message') is None:
        pass
    elif result_query.find('message').text == 'No hits found':
        query_acc = result_query.find('query-title').text
        query_summary = str.splitlines(query_acc) + str.splitlines('No hits found')
        with open(output_file2, "a", encoding="utf-8") as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(query_summary)
        continue
    query_acc = result_query.find('query-title').text
    query_len = result_query.find('query-len').text
    query_info = str.splitlines(query_acc) + str.splitlines(query_len)
    query_hits = result_query.find_all('Hit')
    if mode == 'blastn':
        hit_sciname_total = [] # blastn specific
    else:
        hit_title_total = []   # blastp, blastx
    hit_num_total = []
    for query_hit in query_hits:
        hit_num = query_hit.find('num').text
        hit_num_total.append(hit_num)
        if query_hit.find('accession') is None:
            hit_acc = query_hit.find('id').text
        else:
            hit_acc = query_hit.find('accession').text
        if query_hit.find('title') is None:
            hit_title = "null"
            if mode != 'blastn':
                hit_title_total.append(hit_title)
            else:
                pass
        else:
            hit_title_full = query_hit.find('title').text.split(';') # blastp, blastx
            hit_title = hit_title_full[0:1]                          # blastp, blastx
            if mode != 'blastn':
                hit_title_total = hit_title_total + hit_title        # blastp, blastx
            else:
                pass
        if query_hit.find('taxid') is None:
            hit_taxid = "null"
        else:
            hit_taxid = query_hit.find('taxid').text
        if query_hit.find('sciname') is None:
            hit_sciname = "null"
        else:
            hit_sciname = query_hit.find('sciname').text
        if mode == 'blastn':
            hit_sciname_total.append(hit_sciname) # blastn specific
        else:
            pass
        hit_len = query_hit.find('len').text
        if mode == 'blastn':
            hit_info = str.splitlines(hit_num) + str.splitlines(hit_acc) + str.splitlines(hit_taxid) + str.splitlines(hit_sciname) + str.splitlines(hit_len)             # blastn specific
        elif type(hit_title) == list:
            hit_info = str.splitlines(hit_num) + str.splitlines(hit_acc) + hit_title + str.splitlines(hit_taxid) + str.splitlines(hit_sciname) + str.splitlines(hit_len) # blastp, blastx
        elif type(hit_title) == str:
            hit_info = str.splitlines(hit_num) + str.splitlines(hit_acc) + str.splitlines(hit_title) + str.splitlines(hit_taxid) + str.splitlines(hit_sciname) + str.splitlines(hit_len) # blastp, blastn
        hit_hsps = query_hit.find_all('hsps')
        for hit_hsp in hit_hsps:
            alinum = []
            bitscore = []
            score = []
            evalue = []
            ident = []
            query_f = []
            query_t = []
            subj_f = []
            subj_t = []
            alilen = []
            gaps = []
            for hsp_num in hit_hsp.find_all('num'):
                alinum.append(hsp_num.text)
            for hsp_bitscore in hit_hsp.find_all('bit-score'):
                bitscore.append(hsp_bitscore.text)
            for hsp_score in hit_hsp.find_all('score'):
                score.append(hsp_score.text)
            for hsp_eval in hit_hsp.find_all('evalue'):
                evalue.append(hsp_eval.text)
            for hsp_ident in hit_hsp.find_all('identity'):
                ident.append(hsp_ident.text)
            for hsp_query_f in hit_hsp.find_all('query-from'):
                query_f.append(hsp_query_f.text)
            for hsp_query_t in hit_hsp.find_all('query-to'):
                query_t.append(hsp_query_t.text)
            for hsp_subj_f in hit_hsp.find_all('hit-from'):
                subj_f.append(hsp_subj_f.text)
            for hsp_subj_t in hit_hsp.find_all('hit-to'):
                subj_t.append(hsp_subj_t.text)
            for hsp_alilen in hit_hsp.find_all('align-len'):
                alilen.append(hsp_alilen.text)
            for hsp_gaps in hit_hsp.find_all('gaps'):
                gaps.append(hsp_gaps.text)
            if mode == 'blastn':
                query_s = [] # blastn specific
                subj_s = []  # blastn specific
                for hsp_query_s in hit_hsp.find_all('query-strand'): # blastn specific
                    query_s.append(hsp_query_s.text)                 # blastn specific
                for hsp_subj_s in hit_hsp.find_all('hit-strand'):    # blastn specific
                    subj_s.append(hsp_subj_s.text)                   # blastn specific
            elif mode == 'blastp':
                positive = [] # blastp, blastx
                for hsp_positive in hit_hsp.find_all('positive'):    # blastp, blastx
                    positive.append(hsp_positive.text)               # blastp, blastx
            elif mode == 'blastx':
                positive = [] # blastp, blastx
                query_fr = [] # blastx specific
                for hsp_positive in hit_hsp.find_all('positive'):    # blastp, blastx
                    positive.append(hsp_positive.text)               # blastp, blastx
                for hsp_query_fr in hit_hsp.find_all('query-frame'): # blastx specific
                    query_fr.append(hsp_query_fr.text)               # blastx specific
            else:
                pass
        cont_num = len(list(alinum))
        identity = [x / y * 100 for (x, y) in zip(map(int, ident), map(int, alilen))]
        for i in range(cont_num):
            if mode == 'blastx':
                aa_conv_len = int(query_len) // 3                                # blastx specific
                coverage = ((int(alilen[i]) - int(gaps[i])) / aa_conv_len * 100) # blastx specific
            elif mode != 'blastx':
                coverage = ((int(alilen[i]) - int(gaps[i])) / int(query_len) * 100)
            coverage_f = "{:.0f}".format(coverage)
            identity_f = "{:.0f}".format(identity[i])
            if mode == 'blastn' and type(hit_title) == list:
                tsv_out = query_info + hit_info + alinum[i:i+1] + bitscore[i:i+1] + score[i:i+1] + evalue[i:i+1] + ident[i:i+1] + query_s[i:i+1] + subj_s[i:i+1] + query_f[i:i+1] + query_t[i:i+1] + subj_f[i:i+1] + subj_t[i:i+1] + alilen[i:i+1] + gaps[i:i+1] + str.splitlines(coverage_f) + str.splitlines(identity_f) + hit_title
            elif mode == 'blastn' and type(hit_title) == str:
                tsv_out = query_info + hit_info + alinum[i:i+1] + bitscore[i:i+1] + score[i:i+1] + evalue[i:i+1] + ident[i:i+1] + query_s[i:i+1] + subj_s[i:i+1] + query_f[i:i+1] + query_t[i:i+1] + subj_f[i:i+1] + subj_t[i:i+1] + alilen[i:i+1] + gaps[i:i+1] + str.splitlines(coverage_f) + str.splitlines(identity_f) + str.splitlines(hit_title)
            elif mode == 'blastp':
                tsv_out = query_info + hit_info + alinum[i:i+1] + bitscore[i:i+1] + score[i:i+1] + evalue[i:i+1] + ident[i:i+1] + positive[i:i+1] + query_f[i:i+1] + query_t[i:i+1] + subj_f[i:i+1] + subj_t[i:i+1] + alilen[i:i+1] + gaps[i:i+1] + str.splitlines(coverage_f) + str.splitlines(identity_f)
            elif mode == 'blastx':
                tsv_out = query_info + hit_info + alinum[i:i+1] + bitscore[i:i+1] + score[i:i+1] + evalue[i:i+1] + ident[i:i+1] + positive[i:i+1] + query_fr[i:i+1] + query_f[i:i+1] + query_t[i:i+1] + subj_f[i:i+1] + subj_t[i:i+1] + alilen[i:i+1] + gaps[i:i+1] + str.splitlines(coverage_f) + str.splitlines(identity_f)
            with open(output_file1, "a", encoding="utf-8") as f:
                writer = csv.writer(f, delimiter='\t')
                writer.writerow(tsv_out)
    query_summary = str.splitlines(query_acc) + str.splitlines(hit_num_total[-1])
    if mode == 'blastn':
        total = collections.Counter(hit_sciname_total) # blastn specific
    else:
        total = collections.Counter(hit_title_total)   # blastp, blastx
    with open(output_file2, "a", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(query_summary)
        for k, v in total.items():
            hit_summary = str.splitlines(str(k)) + str.splitlines(str(v))
            writer.writerow(hit_summary)  

print("\nDone\n")
