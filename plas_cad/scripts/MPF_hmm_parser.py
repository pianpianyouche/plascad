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
###################################### parsing MPF module #############################################
def MPF_parsing(in_domtblout, coverage):
    dic = OrderedDict()
    MPF_restult = open(str(in_domtblout) + "_" + "parsed_out", "w")
    wanted1 = set()
    wanted2 = set()
    for line in open(in_domtblout, "r"):
        if not line.startswith("#"):
            KEY = "\t".join(str(line).strip().split()[0:6])
            c_value = str(line).strip().split()[11]
            domain_start = str(line).strip().split()[17]
            domain_end = str(line).strip().split()[18]
            wanted1.add(KEY + "\t" + c_value + "\t" + domain_start + "\t" + domain_end)
    for line in wanted1:
        try:
            items = re.split("\t", line.strip())
            key = '\t'.join(items[:6])
            value1 = items[6]
            value2 = items[7]
            value3 = items[8]
            if key not in dic:
                dic[key] = [[item] for item in items[6:]]
            else:
                if float(value1) < float(dic[key][0][0]):
                    dic[key][0] = value1
                if float(value2) < float(dic[key][1][0]):
                    dic[key][1][0] = value2
                if float(value3) > float(dic[key][2][0]):
                    dic[key][2][0] = value3
        except:
            pass
    for k, v in dic.items():
        wanted2.add('{}\t{}\t{}\t{}\n'.format(k, *map(''.join, (v))))
###################################### parsing based on coverage ############################
    MPF_wanted = OrderedDict()
    for line in wanted2:
        id = str(line).strip().split("\t")[0]
        lis = str(line).strip().split("\t")
        if ((float(lis[8]) - float(lis[7])) / float(lis[5])) * 100 >= int(coverage) or \
                                            ((float(lis[8]) - float(lis[7])) / float(lis[2])) * 100 >= int(coverage) and float(lis[6]) <= 0.01:
            MPF_wanted[id] = line
    for v in MPF_wanted.values():
        MPF_restult.write(v)
    MPF_restult.close()
##################################### MPF classification ####################################
for root, folders, files in os.walk(os.getcwd()):
    for i in fnmatch.filter(files, '*hmm_domtblout'):
        MPF_parsing(i, 50)
##################################### ATPase and T4CP summary ########################################
catATPase = 'cat ' + '*ATPase.hmm_domtblout_parsed_out> ' + str(args.i) + '_ATPase_temp_out'
catT4CP = 'cat ' + '*T4CP*_domtblout_parsed_out> ' + str(args.i) + '_T4CP_temp_out'
os.system(catATPase)
os.system(catT4CP)
#################################### ATPase and T4CP parsing result #######################################
def MPFF_parsed_result(summary, faa):
    result = open(str(summary).rsplit("_", 1)[0] + "_parsed_result_out", "w")
    diclocation = {}
    location_wanted = set()
    for record in SeqIO.parse(open(faa, "r"), "fasta"):
        diclocation[record.id] = record.description.split("#")[1].strip() + "-" + record.description.split("#")[2].strip() + "\t" + record.description.split("#")[3].strip()
    dicATPase_T4CP = {}
    for line in open(summary, "r"):
        key = str(line).strip().split("\t")[0]
        e_value = float(str(line).strip().split("\t")[6])
        if key not in dicATPase_T4CP.keys():
            dicATPase_T4CP[key] = str(line).strip().split("\t")[3] + "\t" + \
                          str(line).strip().split("\t")[6] + "\t" + str(diclocation[key]) + "\n"
        else:
            if float(e_value) < float(str(dicATPase_T4CP[key]).split("\t")[1]):
                dicATPase_T4CP[key] = str(line).strip().split("\t")[3] + "\t" + \
                              str(line).strip().split("\t")[6] + "\t" + str(diclocation[key]) + "\n"
    for key, value in dicATPase_T4CP.items():
        location_wanted.add(str(key) + "\t" + str(value).strip())
    dic_wanted = {}
    for line in location_wanted:
        key = str(line).split("\t")[0].rsplit("_", 1)[0]
        e_value = float(str(line).strip().split("\t")[2])
        if key not in dic_wanted.keys():
            dic_wanted[key] = "\t".join(str(line).split("\t")[1:])
        else:
            if float(e_value) < float(str(dic_wanted[key]).split("\t")[1]):
                dic_wanted[key] = "\t".join(str(line).split("\t")[1:])
    for k, v in dic_wanted.items():
        result.write(str(k) + "\t" + str(v).strip() + "\n")
    result.close()
MPFF_parsed_result(str(args.i) + '_ATPase_temp_out', str(args.i) + "_MOB_temp_mob.faa")
MPFF_parsed_result(str(args.i) + '_T4CP_temp_out', str(args.i) + "_MOB_temp_mob.faa")
################################################################################################
for root, folders, files in os.walk(os.getcwd()):
    for i in fnmatch.filter(files, '*hmm_domtblout_parsed_out'):
        MPFF_parsed_result(i, str(args.i) + "_MOB_temp_mob.faa")
################################################################################################
