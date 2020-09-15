#! /usr/bin/env python3.6
from collections import OrderedDict
import re
from Bio import SeqIO
import os
import fnmatch
import argparse
###################################### Arguments and declarations ########################################
parser = argparse.ArgumentParser()
parser.add_argument("-i",
                    help="input plasmids file prefix for classification",type=str)
parser.add_argument("-db",
                    help="ARGs database", type=str)
parser.add_argument("-db_structure",
                    help="ARGs database structure", type=str)
args = parser.parse_args()
def plasmids_ARGs_filter(blast_out, ARGsdb):
    ARGs_filter = []
    a = {}
    for record in SeqIO.parse(open(ARGsdb, "r"), "fasta"):
        a[str(record.id).strip()] = len(str(record.seq))
    for line in open(blast_out, "r"):
        lis = line.split("\t")
        if float(lis[2]) >= 80 and float(lis[3]) * 100 / float(a[lis[1]]) >= 70:
            ARGs_filter.append(line)
    return ARGs_filter
###################################### description ########################################
def plasmids_des(des):
    ARGs_des = []
    b = {}
    for line in open(des, "r"):
        b[str(line).split("\t")[0]] = str(line).split("\t")[1]
    ARGs = plasmids_ARGs_filter(str(args.i) + "_ARGs_blasp_result_out", str(args.db))
    for i in ARGs:
        if str(i).split("\t")[1] in b.keys():
            ARGs_des.append(str(i).split("\t")[0] + "\t" + b[str(i).split("\t")[1]] + "\t" + str(i).split("\t")[10])
    return ARGs_des
# # ###################################### ARGs parsing result ########################################
def plasmids_ARGs_parser(faa):
    diclocation = {}# location and oritation
    result = open(str(args.i) + "_ARGs_annotation_out", "w")
    for record in SeqIO.parse(open(faa, "r"), "fasta"):
        diclocation[record.id] = record.description.split("#")[1].strip() + "-" + record.description.split("#")[2].strip() + "\t" + record.description.split("#")[3].strip()
    plas_des = plasmids_des(str(args.db_structure))
    for line in plas_des:
        search = str(line).split("\t")[0]
        if search in diclocation.keys():
            result.write(str(search).rsplit("_",1)[0] + "\t" + str(line).split("\t")[1].strip() + "\t" + str(line).split("\t")[2].strip() + "\t" + str(diclocation[search]).strip() + "\n")
    result.close()

if __name__ == "__main__":
    plasmids_ARGs_parser(str(args.i) + ".faa")










