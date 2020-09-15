#! /usr/bin/env python3.6
from collections import OrderedDict
import shutil
import os
import glob
import argparse
import sys
import subprocess
sys.path.append("..")
import plas_cad

def main():
    usage = ("usage: plas_cad -i your.plasmid.seqs.fasta")
    version = 'plas_cad {v}'.format(v=plas_cad.__version__)
###################################### checking dependencies ########################################
    list_cmd = ['prodigal', 'blastp', 'hmmsearch']
    for cmd in list_cmd:
        exist = subprocess.call('command -v '+ cmd + '>> /dev/null', shell=True)
        if exist == 0:
            pass
        else:
            print(cmd + " not exist in path!")
            sys.exit()
###################################### Arguments and declarations ########################################
    parser = argparse.ArgumentParser()
    parser.add_argument("-i",
                        help="input plasmids file for classification", type=str,
                        default='example/example.fasta')
    parser.add_argument("-cMOBB",
                        help="alignment coverage for MOBB HMM profile",
                        default=75)
    parser.add_argument("-cMOBC",
                        help="alignment coverage for MOBC HMM profile",
                        default=75)
    parser.add_argument("-cMOBF",
                        help="alignment coverage for MOBF HMM profile",
                        default=75)
    parser.add_argument("-cMOBT",
                        help="alignment coverage for MOBT HMM profile",
                        default=75)
    parser.add_argument("-cMOBPB",
                        help="alignment coverage for MOBPB HMM profile",
                        default=75)
    parser.add_argument("-cMOBH",
                        help="alignment coverage for MOBH HMM profile",
                        default=70)
    parser.add_argument("-cMOBP",
                        help="alignment coverage for MOBP HMM profile",
                        default=65)
    parser.add_argument("-cMOBV",
                        help="alignment coverage for MOBV HMM profile",
                        default=60)
    parser.add_argument("-cMOBQ",
                        help="alignment coverage for MOBQ HMM profile",
                        default=55)
    args = parser.parse_args()
    directory =  os.path.abspath(os.path.dirname(__file__))
    file_name, file_ext = os.path.splitext(args.i)
###################################### Prodigal ###########################################################
    cmdprodigal = "prodigal"  + " -i "  + str(args.i) + \
                " -a " + str(file_name) + ".faa " + " -q -o temp.txt"

    os.system(cmdprodigal)
    os.remove("temp.txt")
###################################### MOB hmmer ###############################################################
    mob_hmm = os.path.join(directory, "database/hmm_module/MOB_hmm")
    for root, dirnames, filenames in os.walk(mob_hmm):
        for i in filenames:
            cmdhmmsearch = "hmmsearch" + " --domtblout " + str(file_name) + "_" + str(i) + \
                       "_domtblout " + str(os.path.join(mob_hmm, i)) + " " + str(file_name) + ".faa "
            os.system(cmdhmmsearch)
###################################### hmmer parsing #########################################################
    cmdmobparsing = 'python ' + os.path.join(directory, "scripts/MOB_parser.py") + ' -i ' + str(file_name) + ' -cMOBB ' + str(args.cMOBB) + ' -cMOBC ' + str(args.cMOBC) + \
        ' -cMOBF ' + str(args.cMOBF)+ ' -cMOBT ' + str(args.cMOBT) + ' -cMOBPB ' + str(args.cMOBPB) + ' -cMOBH ' + str(args.cMOBH) \
        + ' -cMOBP ' + str(args.cMOBP) + ' -cMOBV ' + str(args.cMOBV) + ' -cMOBQ ' + str(args.cMOBQ) + "\n"
    os.system (cmdmobparsing)
###################################### MPF system hmmer #######################################################
    MPF_hmm = os.path.join(directory, "database/hmm_module/MPF_system_hmm")
    for root, dirnames, filenames in os.walk(MPF_hmm):
        for i in filenames:
            cmdMPF = "hmmsearch" + " --domtblout " + str(file_name) + "_mob_" + str(i) + \
                "_domtblout " + str(os.path.join(MPF_hmm, i)) + " " + str(file_name) + "_MOB_temp_mob.faa"
            os.system(cmdMPF)
###################################### MPF parsing ##########################################################
    cmdATPase_T4CP_parsing = 'python ' + os.path.join(directory, 'scripts/MPF_hmm_parser.py') + ' -i ' + str(file_name) + "\n"
    os.system(cmdATPase_T4CP_parsing)
###################################### MPFF summary ##########################################################
    cmdcatmpff = 'cat ' + str(file_name) + '_ATPase_temp_parsed_result_out ' + str(file_name) +  '_T4CP_temp_parsed_result_out ' \
             + str(file_name) + "*MPFF*_result_out > " + str(file_name) + "_MPFF_summary_out"
    cmdcatmpfg = 'cat ' + str(file_name) + '_ATPase_temp_parsed_result_out ' + str(file_name) +  '_T4CP_temp_parsed_result_out ' \
             + str(file_name) + "*MPFG*_result_out > " + str(file_name) + "_MPFG_summary_out"
    cmdcatmpfi = 'cat ' + str(file_name) + '_ATPase_temp_parsed_result_out ' + str(file_name) +  '_T4CP_temp_parsed_result_out ' \
             + str(file_name) + "*MPFI*_result_out > " + str(file_name) + "_MPFI_summary_out"
    cmdcatmpft = 'cat ' + str(file_name) + '_ATPase_temp_parsed_result_out ' + str(file_name) +  '_T4CP_temp_parsed_result_out ' \
             + str(file_name) + "*MPFT*_result_out > " + str(file_name) + "_MPFT_summary_out"
    os.system(cmdcatmpff)
    os.system(cmdcatmpfg)
    os.system(cmdcatmpfi)
    os.system(cmdcatmpft)
###################################### MPFF classification ##########################################################
    cmdConjparsing = 'python ' + os.path.join(directory, 'scripts/Conj_parser.py') + ' -i ' + str(file_name) + "\n"
    os.system(cmdConjparsing)
    cmdcatconj = 'cat ' + str(file_name) + "*Conj_summary_out > " + str(file_name) + "_Conj_out"
    os.system(cmdcatconj)
###################################### Plasmids classification ######################################################
    cmdPlasparsing = 'python ' + os.path.join(directory, 'scripts/Conj_mob_classification.py') + ' -i ' + str(file_name) + "\n"
    os.system(cmdPlasparsing)
###################################### Plasmids ARGs identification ######################################################
    cmdblastp = "blastp" + " -query " + str(file_name) + ".faa" + " -db "\
    + os.path.join(directory, 'database/ARGsdb/ARGsDB') + " -outfmt 6 -evalue 1e-5 -num_threads 20 " +\
            "| sort -k1,1 -k12,12nr -k11,11n | sort -u -k1,1 --merge > " + str(file_name) + "_ARGs_blasp_result_out"
    os.system(cmdblastp)
###################################### Plasmids ARGs parsing ######################################################
    cmdARGsparsing = 'python ' + os.path.join(directory, 'scripts/plasmids_ARGs_parser.py') + ' -i ' + str(file_name) + \
    ' -db ' + os.path.join(directory, 'database/ARGsdb/ARGsDB.fasta') + ' -db_structure ' + os.path.join(directory, 'database/ARGsdb/ARGsDB_des.txt')
    os.system(cmdARGsparsing)
###################################### Plasmids ARGs summary ######################################################
    cmdARGssummary = 'python ' + os.path.join(directory, 'scripts/plasmids_result_summary.py') + ' -i ' + str(file_name)
    os.system(cmdARGssummary)
###################################### summary result ######################################################
    cmdcatsummary = 'cat ' + str(file_name) + '_Conj_plasmids_ids_result_out ' + str(file_name) + '_mob_unconj_plasmids_ids_result_out ' +\
    str(file_name) + '_unmob_plasmids_ids_result_out > ' + str(file_name) + "_plasmids_classification_sum.txt"
    os.system(cmdcatsummary)
###################################### plasmids map plotting ######################################################
    cmdplot = 'python ' + os.path.join(directory, 'scripts/Plasmids_plot.py') + ' -i ' +  str(file_name)
    os.system(cmdplot)
###################################### remove temp  ######################################################
    os.system('rm -rf *out')
    os.system('rm -rf *faa')

###################################### function  ######################################################
if __name__ == '__main__':
    main()
