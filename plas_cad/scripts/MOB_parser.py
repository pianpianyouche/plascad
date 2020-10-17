#! /usr/bin/env python3.6
from collections import OrderedDict
import re
from Bio import SeqIO
import os
import argparse
###################################### Arguments and declarations ########################################
parser = argparse.ArgumentParser()
parser.add_argument("-i",
                    help="input plasmids file prefix for classification",type=str)
parser.add_argument("-cMOBB",
                    help="alignment coverage for MOBB HMM profile",
                    type=int)
parser.add_argument("-cMOBC",
                    help="alignment coverage for MOBC HMM profile",
                    type=int)
parser.add_argument("-cMOBF",
                    help="alignment coverage for MOBF HMM profile",
                    type=int)
parser.add_argument("-cMOBT",
                    help="alignment coverage for MOBT HMM profile",
                    type=int)
parser.add_argument("-cMOBPB",
                    help="alignment coverage for MOBPB HMM profile",
                    type=int)
parser.add_argument("-cMOBH",
                    help="alignment coverage for MOBH HMM profile",
                    type=int)
parser.add_argument("-cMOBP",
                    help="alignment coverage for MOBP HMM profile",
                    type=int)
parser.add_argument("-cMOBV",
                    help="alignment coverage for MOBV HMM profile",
                    type=int)
parser.add_argument("-cMOBQ",
                    help="alignment coverage for MOBQ HMM profile",
                    type=int)
args = parser.parse_args()
#print ("MOB_parser.py was developed by You Che")
###################################### parsing each MOB module #############################################
def MOB_parsing(in_domtblout, coverage):
    dic = OrderedDict()
    MOB_restult = open(str(in_domtblout) + "_" + "parsed_out", "w")
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
#################################### parsing based on coverage ############################
    mob_wanted = OrderedDict()
    for line in wanted2:
        id = str(line).strip().split("\t")[0]
        lis = str(line).strip().split("\t")
        if ((float(lis[8]) - float(lis[7])) / float(lis[5])) * 100 >= int(coverage) or \
                                        ((float(lis[8]) - float(lis[7])) / float(lis[2])) * 100 >= int(coverage) and float(lis[6]) <= 0.01:
            mob_wanted[id] = line
    for v in mob_wanted.values():
        MOB_restult.write(v)
    MOB_restult.close()
#################################### MOB classification ####################################
MOB_parsing(str(args.i) + "_MOBB.hmm_domtblout", args.cMOBB)
MOB_parsing(str(args.i) + "_MOBC.hmm_domtblout", args.cMOBC)
MOB_parsing(str(args.i) + "_MOBF.hmm_domtblout", args.cMOBF)
MOB_parsing(str(args.i) + "_MOBT.hmm_domtblout", args.cMOBT)
MOB_parsing(str(args.i) + "_MOBPB.hmm_domtblout", args.cMOBPB)
MOB_parsing(str(args.i) + "_MOBH.hmm_domtblout", args.cMOBH)
MOB_parsing(str(args.i) + "_MOBP.hmm_domtblout", args.cMOBP)
MOB_parsing(str(args.i) + "_MOBV.hmm_domtblout", args.cMOBV)
MOB_parsing(str(args.i) + "_MOBQ.hmm_domtblout", args.cMOBQ)
#################################### summary mob ##############################################
catcmd = 'cat ' + '*MOB*_domtblout_parsed_out > ' + str(args.i) + '_MOB_temp_out'
os.system(catcmd)
#################################### mob parsing result #######################################
def mob_classification(summary_mob, faa):
    location_wanted = set()
    diclocation = {}
    for record in SeqIO.parse(open(faa, "r"), "fasta"):
        diclocation[record.id] = record.description.split("#")[1].strip() + "-" + record.description.split("#")[2].strip() + "\t" + record.description.split("#")[3].strip()
    dicmob = {}
    for line in open(summary_mob, "r"):
        key = str(line).strip().split("\t")[0]
        e_value = float(str(line).strip().split("\t")[6])
        if key not in dicmob.keys():
            dicmob[key] = str(line).strip().split("\t")[3] + "\t" +\
                          str(line).strip().split("\t")[6] + "\t" + str(diclocation[key]) + "\n"
        else:
            if float(e_value) < float(str(dicmob[key]).split("\t")[1]):
                dicmob[key] = str(line).strip().split("\t")[3] + "\t" +\
                          str(line).strip().split("\t")[6] + "\t" + str(diclocation[key]) + "\n"
    for k, v in dicmob.items():
        location_wanted.add('{}\t{}\n'.format(str(k).strip(), str(v).strip()))
    dic_wanted = {}
    result = open(str(summary_mob).rsplit("_", 1)[0] + "_parsed_result_out", "w")
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
    mob_result = open(str(summary_mob).rsplit("_", 1)[0] + "_mob.faa" , "w")
    unmob_result = open(str(summary_mob).rsplit("_", 1)[0] + "_unmob.faa" , "w")
    for record in SeqIO.parse(open(faa, "r"), "fasta"):
        if str(record.id).rsplit("_", 1)[0] in dic_wanted.keys():
            SeqIO.write(record, mob_result, "fasta")
        else:
            SeqIO.write(record, unmob_result, "fasta")
    result.close()
    mob_result.close()
    unmob_result.close()
if __name__ == "__main__":
    mob_classification(str(args.i) + '_MOB_temp_out', str(args.i) + ".faa")
#################################### mob unmob faa extraction #######################################
