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
def result_summary(ARGs, loc, ids, p_type):
    result_loc = open(str(loc).rsplit("_",2)[0] + "_loc_sum.txt" , "w")
    result_ids = open(str(ids).rsplit("_",2)[0] + "_ids_result_out", "w")
    id_lis = {}
    loc_lis = []
    for i in open(loc, "r"):
        loc_lis.append(str(i).strip())
    for i in open(ids, "r"):
        id_lis[str(i).strip()] = []
    for i in open(ARGs, "r"):
        search = str(i).strip().split("\t")[0]
        if search in id_lis.keys():
            id_lis[search].append(str(i).split("\t")[1])
            loc_lis.append(str(i).strip())
    loc_lis_sorted = sorted(loc_lis, key=lambda x: (str(x).strip().split()[0]))
    for i in loc_lis_sorted:
        result_loc.write(str(i) + "\n")
    for k, v in id_lis.items():
        result_ids.write(str(k).strip() + "\t" + p_type + "\t" + "|".join(v)+ "\n")
    result_loc.close()
    result_ids.close()
###################################### unmob plasmids ########################################
def result_summary_unmob(ARGs, ids, p_type):
    result_ids = open(str(ids).rsplit("_",2)[0] + "_ids_result_out", "w")
    id_lis = {}
    for i in open(ids, "r"):
        id_lis[str(i).strip()] = []
    for i in open(ARGs, "r"):
        search = str(i).strip().split("\t")[0]
        if search in id_lis.keys():
            id_lis[search].append(str(i).split("\t")[1])
    for k, v in id_lis.items():
        result_ids.write(str(k).strip() + "\t" + p_type + "\t" + "|".join(v)+ "\n")
    result_ids.close()
###################################### unmob plasmids ########################################
if __name__ == "__main__":
    result_summary(str(args.i) + "_ARGs_annotation_out", str(args.i) + "_Conj_plasmids_loc_out", str(args.i) + "_Conj_plasmids_id_out", "Conj")
    result_summary(str(args.i) + "_ARGs_annotation_out", str(args.i) + "_mob_unconj_plasmids_loc_out", str(args.i) + "_mob_unconj_plasmids_id_out", "mob_unconj")
    result_summary_unmob(str(args.i) + "_ARGs_annotation_out", str(args.i) + "_unmob_plasmids_id_out", "unmob")










