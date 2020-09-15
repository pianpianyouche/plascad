#! /usr/bin/env python3
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
def Conj_mob_classification(conj_in, mob_in):
    Conj_sum = open(str(conj_in).rsplit("_",2)[0] + "_Conj_plasmids_loc_out", "w")
    Conj_id = open(str(conj_in).rsplit("_",2)[0] + "_Conj_plasmids_id_out", "w")
    mob_unconj_id = open(str(conj_in).rsplit("_",2)[0] + "_mob_unconj_plasmids_id_out", "w")
    mob_unconj = open(str(conj_in).rsplit("_",2)[0] + "_mob_unconj_plasmids_loc_out", "w")
    Conjugative = []
    Mob_unconj = set()
    Conjugative_id = set()
    for line in open(conj_in, "r"):
        Conjugative.append(line)
    for line in Conjugative:
        Conjugative_id.add(str(line).split("\t")[0].strip())
    for line in open(mob_in, "r"):
        if str(line).split("\t")[0] in Conjugative_id:
            Conjugative.append(line)
        else:
            Mob_unconj.add(line)
    Conjugative_sorted = sorted(Conjugative, key=lambda x: (str(x).strip().split()[0]))
    for i in Conjugative_sorted:
        Conj_sum.write(i)
    for i in Conjugative_id:
        Conj_id.write(i + "\n")
    for i in Mob_unconj:
        mob_unconj.write(i)
        mob_unconj_id.write(str(i).strip().split("\t")[0] + "\n")
    Conj_sum.close()
    Conj_id.close()
    mob_unconj.close()
    mob_unconj_id.close()
#Conj_mob_classification(str(args.i) + "_Conj_out", str(args.i) + "_MOB_temp_parsed_result_out")
###################################### unmob id extraction ########################################
def unmob_id(mob_result, all_fasta):
    result = open(str(args.i) + "_unmob_plasmids_id_out", "w")
    all_mob_id = set()
    for line in open(mob_result, "r"):
        all_mob_id.add(str(line).split("\t")[0])
    with open(all_fasta, "r") as fdb:
        for record in SeqIO.parse(fdb, "fasta"):
            if str(record.id) not in all_mob_id:
                result.write(str(record.id).strip() + "\n")
    result.close()
###################################### running  ########################################
if __name__ == "__main__":
    Conj_mob_classification(str(args.i) + "_Conj_out", str(args.i) + "_MOB_temp_parsed_result_out")    
    unmob_id(str(args.i) + "_MOB_temp_parsed_result_out", str(args.i) + ".fasta")
#print "Plasmids classification finished!!!"
