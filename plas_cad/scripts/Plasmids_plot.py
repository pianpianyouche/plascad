#! /usr/bin/env python3.6
from xml.etree.ElementTree import Element, SubElement, Comment, tostring
from ElementTree_pretty import prettify
from Bio import SeqIO
import argparse
import os
###################################### Arguments and declarations ########################################
parser = argparse.ArgumentParser()
parser.add_argument("-i",
                    help="input plasmids file prefix for classification",type=str)
args = parser.parse_args()

####################################################### plasmids plot process #####################################
curdir = os.getcwd()
directory =  os.path.abspath(os.path.dirname(__file__))
result_dir = os.path.join(curdir, str(args.i) +'_Conjugative_plasmids_map')
if not os.path.exists(result_dir):
    os.makedirs(result_dir)
cmdcpjs = 'cp -r ' + os.path.join(directory, 'js ') + result_dir
os.system (cmdcpjs)
####################################################### main_program #####################################
def html_prepare(lis_in):
    m = 33
    n = 33
    length = lis_in[0].split("\t")[5]
    ID = lis_in[0].split("\t")[0]
    tick = str(round(float(length)/50))
    tick2 = str(round(float(length)/12))
    result = open(os.path.join(result_dir, ID + ".html"), "w")
    root = Element("html")
    root.set('version', '1.0')
    head = SubElement(root, 'head')
    title = SubElement(head, 'script', attrib= {"src":"js/angularplasmid.complete.min.js"})
    title.text = " "
# title.text = " src='js/angularplasmid.complete.min.js'"
    body = SubElement(root, 'body')
    style = SubElement(body, 'style')
    ti_style = SubElement(style, "setting")
    ti_style.text = " body {font-family: 'Lato';font-weight:400;}" \
                ".boundary {stroke-dasharray:2,2;stroke-width:2px}" \
                ".mdlabel {font-size:14px}" \
                ".smlabel {font-size:8px}" \
                ".white {fill:#fff}" \
                ".red {fill:rgb(192,64,64)}" \
                ".purple {fill:rgb(192,64,192)}" \
                ".blue {fill:rgb(64,192,192)}" \
                ".green {fill:rgb(64,192,64)}" \
                ".labelline {stroke:#333;stroke-dasharray:2,2;stroke-width:2px;}" \
                ".gold {fill:rgb(192,128,64)}" \
                ""
    plasmid = SubElement(body, 'plasmid', attrib={"sequencelength": length, "plasmidheight": '700', "plasmidwidth":'700'})
    plasmidtrack =  SubElement(plasmid, 'plasmidtrack', attrib={"trackstyle":"fill:#ccc", "width":"5", "radius":"150"})
    plasmidtrack.text = " "
    plasmidtrack2 =  SubElement(plasmid, 'plasmidtrack', attrib={"trackstyle":"fill:rgba(225,225,225,0.5)","radius":"140"})
    tracklabel =  SubElement(plasmidtrack2, "tracklabel", attrib={"text":ID, "labelstyle":"font-size:20px;font-weight:400"})
    tracklabel.text = " "
    tracklabel =  SubElement(plasmidtrack2, "tracklabel", attrib={"text":length + " bp", "labelstyle":"ffont-size:10px", "vadjust":"20"})
    tracklabel.text = " "
    trackscale = SubElement(plasmidtrack2, "trackscale", attrib={"interval":tick, "style":"stroke:#999", "ticksize":"3"})
    trackscale.text = " "
    trackscale = SubElement(plasmidtrack2, "trackscale", attrib={"interval":tick, "style":"stroke:#999","direction": "in", "ticksize":"3"})
    trackscale.text = " "
    trackscale = SubElement(plasmidtrack2, "trackscale", attrib={"interval":tick2, "style":"stroke:#f00", "direction":"in", "showlabels":"1", "labelstyle":"fill:#999;stroke:none;text-anchor:middle;alignment-baseline:middle;font-size:10px"})
    trackscale.text = " "

    for i in lis_in:
        if "MOB" in i:
            mob_type = str(i).split("\t")[1]
            start = str(i).split("\t")[3].split("-")[0]
            end = str(i).split("\t")[3].split("-")[1]
            if str(i).split("\t")[4] == str(1):
                arrow_s = str(-2)
                arrow_e = str(2)
            else:
                arrow_s = str(2)
                arrow_e = str(-2)
            trackmarker = SubElement(plasmidtrack2, "trackmarker", attrib={"start":start, "end":end,"markerstyle":"fill:rgba(85,0,170,0.9)", "arrowendlength":arrow_e, "arrowstartlength":arrow_s})
            markerlabel = SubElement(trackmarker, "markerlabel", attrib={"type":"path", "class":"mdlabel purple","valign":"outer","vadjust":"23","text":mob_type})
            markerlabel.text = " "
            trackmarker = SubElement(plasmidtrack2, "trackmarker", attrib={"start":start, "end":end,"markerstyle":"fill:rgba(238,221,255,0.6)", "wadjust":"-5", "vadjust":"25"})
            trackmarker.text = " "

        if "ATPase" in i:
            start = str(i).split("\t")[3].split("-")[0]
            end = str(i).split("\t")[3].split("-")[1]
            if str(i).split("\t")[4] == str(1):
                arrow_s = str(-2)
                arrow_e = str(2)
            else:
                arrow_s = str(2)
                arrow_e = str(-2)
            trackmarker = SubElement(plasmidtrack2, "trackmarker", attrib={"start":start, "end":end,"markerstyle":"fill:rgba(0,85,170,0.9)", "arrowendlength":arrow_e, "arrowstartlength":arrow_s})
            markerlabel = SubElement(trackmarker, "markerlabel", attrib={"type":"path", "class":"mdlabel blue","valign":"outer","vadjust":"23","text":"ATPase"})
            markerlabel.text = " "
            trackmarker = SubElement(plasmidtrack2, "trackmarker", attrib={"start":start, "end":end,"markerstyle":"fill:rgba(221,238,255,0.6)", "wadjust":"-5", "vadjust":"25"})
            trackmarker.text = " "

        if "T4CP" in i:
            start = str(i).split("\t")[3].split("-")[0]
            end = str(i).split("\t")[3].split("-")[1]
            if str(i).split("\t")[4] == str(1):
                arrow_s = str(-2)
                arrow_e = str(2)
            else:
                arrow_s = str(2)
                arrow_e = str(-2)
            trackmarker = SubElement(plasmidtrack2, "trackmarker", attrib={"start":start, "end":end,"markerstyle":"fill:rgba(170,85,0,0.9)", "arrowendlength":arrow_e, "arrowstartlength":arrow_s})
            markerlabel = SubElement(trackmarker, "markerlabel", attrib={"type":"path", "class":"mdlabel gold", "valign":"outer","vadjust":"23", "text":"T4CP"})
            markerlabel.text = " "
            trackmarker = SubElement(plasmidtrack2, "trackmarker", attrib={"start":start, "end":end,"markerstyle":"fill:rgba(255,238,221,0.6)", "wadjust":"-5", "vadjust":"25"})
            trackmarker.text = " "

        if "MPFG" not in i and "MPFI" not in i and "MPFT" not in i and "MPFF" not in i and "MOB" not in i:
            i = i + "\t" + str(n)
            n +=15
            start = str(i).split("\t")[3].split("-")[0]
            end = str(i).split("\t")[3].split("-")[1]
            des = str(i).split("\t")[1].split("__")[1]
            adjARG = str(i).split("\t")[6]
            if str(i).split("\t")[4] == str(1):
                arrow_s = str(-2)
                arrow_e = str(2)
            else:
                arrow_s = str(2)
                arrow_e = str(-2)
            trackmarker = SubElement(plasmidtrack2, "trackmarker", attrib={"start":start, "end":end,"markerstyle":"fill:rgba(170,0,85,0.9)", "arrowendlength":arrow_e, "arrowstartlength":arrow_s})
            markerlabel = SubElement(trackmarker, "markerlabel", attrib={"type":"path", "class":"mdlabel red", "valign":"outer","vadjust":adjARG,"text":des,"showline":"1","lineclass":"labelline"})
            markerlabel.text = " "
            trackmarker = SubElement(plasmidtrack2, "trackmarker", attrib={"start":start, "end":end,"markerstyle":"fill:rgba(255,221,238,0.6)", "wadjust":"-5", "vadjust":"25"})
            trackmarker.text = " "


        if "virB" in i or "traC" in i or "traE" in i or "traH" in i or "traK" in i or "MPFG_41" in i or "MPFG_44" in i or "MPFG_51" in i or "MPFG_52" in i\
                or "traL" in i or "traN" in i or "traU" in i or "traV" in i or "traW" in i or "traI" in i or "traQ" in i or "traM" in i or "traP" in i\
                or "traR" in i or "traY" in i:
            i = i + "\t" + str(m)
            m +=15
            start = str(i).split("\t")[3].split("-")[0]
            end = str(i).split("\t")[3].split("-")[1]
            des = str(i).split("\t")[1].rsplit("_",1)[1]
            adj = str(i).split("\t")[6]
            if str(i).split("\t")[4] == str(1):
                arrow_s = str(-2)
                arrow_e = str(2)
            else:
                arrow_s = str(2)
                arrow_e = str(-2)
            trackmarker = SubElement(plasmidtrack2, "trackmarker", attrib={"start":start, "end":end,"markerstyle":"fill:rgba(85,170,0,0.9)", "arrowendlength":arrow_e, "arrowstartlength":arrow_s})
            markerlabel = SubElement(trackmarker, "markerlabel", attrib={"type":"path", "class":"mdlabel green","valign":"outer","vadjust":adj, "text":des, "showline":"1","lineclass":"labelline"})
            markerlabel.text = " "
            trackmarker = SubElement(plasmidtrack2, "trackmarker",attrib={"start": start, "end": end, "markerstyle": "fill:rgba(238,255,221,0.6)","wadjust": "-5", "vadjust": "25"})
            trackmarker.text = " "
    result.write(prettify(root))
# ####################################################### prepare the ploting plasmids #####################################
def plasmids_prepare(plas_fasta,conj_sum):
    plas_length = {}
    for record in SeqIO.parse(open(plas_fasta, "r"), "fasta"):
        plas_length[str(record.id).strip()] = len(str(record.seq))
    plot = []
    dic = {}
    for line in open(conj_sum, "r"):
        key = str(line).strip().split("\t")[0]
        if key in plas_length:
            line = str(line).strip() + "\t" + str(plas_length[key]).strip()
        if key not in dic.keys():
            dic[key] = [str(line).strip()]
        else:
            dic[key].append(str(line).strip())
    for k, v in dic.items():
        plot.append(v)
    for i in plot:
        html_prepare(i)

if __name__ == "__main__":
    plasmids_prepare(str(args.i) + ".fasta",str(args.i) + "_Conj_plasmids_loc_sum.txt")
