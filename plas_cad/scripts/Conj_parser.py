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
args = parser.parse_args()
def Conjugation_parse(in_summary, num):
    wanted = set()
    ATPase = []
    T4CP = []
    MPFF = []
    MPFF_combined = set()
    Conj_temp = []
    Conj_summary = open(str(in_summary).rsplit("_",1)[0] + "_Conj_summary_out", "w")
    # Conj_id = open(str(in_summary).rsplit("_",1)[0] + "_Conj_id.txt", "w")
    for line in open(in_summary, "r"):
        MPFF_combined.add(line)
        if "_ATPase" in line:
            ATPase.append(str(line).split("\t")[0])
        if "T4CP" in line:
            T4CP.append(str(line).split("\t")[0])
        elif "_ATPase" not in line and "T4CP" not in line:
            MPFF.append(str(line).split("\t")[0])
    for i in MPFF:
        if MPFF.count(i) >= num and i in ATPase and i in T4CP:
            wanted.add(i)
    for line in MPFF_combined:
        if str(line).strip().split("\t")[0] in wanted:
            Conj_temp.append(str(line).strip())
    sorted_Conj_temp = sorted(Conj_temp, key=lambda x: (str(x).strip().split()[0], str(x).strip().split()[1].split("_")[1]))
    for i in sorted_Conj_temp:
        Conj_summary.write(i + "\n")
    # for i in wanted:
    #     Conj_id.write(i + "\n")

if __name__ == "__main__":
    Conjugation_parse(str(args.i) + "_MPFF_summary_out", 5)
    Conjugation_parse(str(args.i) + "_MPFG_summary_out", 4)
    Conjugation_parse(str(args.i) + "_MPFI_summary_out", 4)
    Conjugation_parse(str(args.i) + "_MPFT_summary_out", 3)
