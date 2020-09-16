# Plascad
Plascad is a computationally efficient tool designed for automated plasmid classification, ARGs annotation and plasmid visualization.


# Contents

* [Plascad workflow](#1)
* [Requirements](#2)
* [Installation](#3)
* [Usage](#4)
* [Example](#5)
* [Output files](#6)
* [Tips for visualization](#7)
* [Acknowledgement](#8)
* [Citation](#9)
* [Contact](#10)

<h2 id="1">Plascad workflow</h2>
Plascad first predicts ORFs in the query plasmid sequences and then detects the homologs of genes associated with plasmid transfer by hmmsearch with the built HMM protein profiles. Then, it identifies ARGs by searching against the structed ARGs database (SARGs), at last, a plasmid visualization component using AngularJS is integrated into this pipeline for the visualization of the classified conjugative plasmids.


![](https://github.com/pianpianyouche/plascad/blob/master/Plascad.jpg) 

<h2 id="2">Requirements</h2>

Linux
[Python >=3.6](https://www.python.org/downloads/)  
[biopython](https://biopython.org/)
[Prodigal](https://github.com/hyattpd/Prodigal)
[blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
[hmmer](http://hmmer.org/)


<h2 id="3">Installation</h2>
Plascad can be installed either through conda or pip, though we advise to use Conda. 


Conda

Use [miniconda](https://docs.conda.io/en/latest/miniconda.html) or [anaconda](https://www.anaconda.com/) to install Plascad

conda create -n Plascad -y -c pianpianyouche plascad
conda activate Plascad

pip3

If you have the dependencies (Python >=3.6, blast >=2.7.1, prodigal >=2.6.3, hmmer >=3.2.1) in your PATH, you can install with pip3

pip3 install Plascad

<h2 id="4">Usage</h2>

`usage: Plascad [-h] [-i I] [-cMOBB CMOBB] [-cMOBC CMOBC] [-cMOBF CMOBF]
               [-cMOBT CMOBT] [-cMOBPB CMOBPB] [-cMOBH CMOBH] [-cMOBP CMOBP]
               [-cMOBV CMOBV] [-cMOBQ CMOBQ]`  

Help:  
    `-h, Show this help message and exit`   
    `-i, FASTA file of plasmid sequences`  
    `-cMOB[B,C,F,T,PB,H,P,V,Q], alignment coverage for MOB HMM profile`

<h2 id="5">Example</h2>

curl -OL https://github.com/pianpianyouche/plascad/raw/master/plas_cad/example/example.fasta
Plascad -i example.fasta

<h2 id="6">Output files</h2>

The final output of Plascad includes four files for the given plasmid sequences:  
1) Summary file of the plasmid classification results.  
2) Genetic location of the genes associated with plasmid transfer and antibiotic resistance genes  
      for both mobilizable and conjugative plasmids. 
3) Html files containing the plasmid maps for all the identified conjugative plasmids.

# Output of the example fasta file
1) `example_plasmids_classification_sum.txt` # sumary of the plasmid classification results  

| Name | plasmid type | ARGs |
| :-: | :-: | :-: |
| AJ627386.1 | Conj| chloramphenicol__catA;tetracycline__tetB;tetracycline__tetC;tetracycline__tetD;beta-lactam__TEM-1 |
| NC_002377.1 | Conj |  |
| NC_002483.1 | Conj |  |
| NC_005014.1 | Conj | tetracycline__tetD;tetracycline__tetC;tetracycline__tetB;aminoglycoside__aph(3'')-I;aminoglycoside__aph(6)-I |  

#note: ARGs are displayed as type_subtype  

2) `example_Conj_plasmids_loc_sum.txt` # Genetic location of the genes associated with transfer and the detected ARGs  

| Name | Marker genes | c-value & e-value | Genetic location | strand |
| :-: | :-: | :-: | :-: | :-: |
| AJ627386.1 | MPFG_41 | 1.2e-88 | 35289-35918 | 1 |
| AJ627386.1 | MPFG_44 | 9.9e-53 | 38227-38634 | 1 |
| AJ627386.1 | MPFG_51 | 8.1e-147 | 44998-46998 | -1 |
| AJ627386.1 | MPFG_52 | 9e-123 | 47013-47954 | -1 |
| AJ627386.1 | MPFG_ATPase | 3.8e-178 | 38650-41520 | 1 |
| AJ627386.1 | MPFG_T4CP_1 | 4.4e-193 | 30481-32730 | 1 |
| AJ627386.1 | MOBH | 3.9e-67 | 55497-57404 | 1 |
| AJ627386.1 | chloramphenicol__catA | 1.94e-164 | 17570-18211 | -1 |
| AJ627386.1 | tetracycline__tetB | 0.0 | 22897-24102 | 1 |
| AJ627386.1 | tetracycline__tetC | 2.36e-166 | 24215-24883 | -1 |
| AJ627386.1 | tetracycline__tetD | 1.18e-100 | 24896-25312 | 1 |
| AJ627386.1 | beta-lactam__TEM-1 | 0.0 | 50985-51845 | -1 |
| NC_002377.1 | MPFT_ATPase | 1.2e-276 | 152657-155026 | 1 |
| NC_002377.1 | MPFT_T4CP_1 | 6.9e-155 | 170649-172619 | 1 |
| NC_002377.1 | MPFT_virB3 | 1.6e-29 | 152331-152657 | 1 |
| NC_002377.1 | MPFT_virB6 | 4.7e-55 | 155800-156687 | 1 |
| NC_002377.1 | MPFT_virB8 | 1.7e-76 | 156815-157588 | 1 |
| NC_002377.1 | MPFT_virB9 | 6.3e-97 | 157585-158466 | 1 |
| NC_002377.1 | MOBQ | 5.7e-104 | 30831-34133 | 1 |  

#note: c-value for hmmsearch, e-value for blastp 


3) `example_Conjugative_plasmids_map` # folder containing the maps for all the identified conjugative plasmids  

![](https://github.com/pianpianyouche/plascad/blob/master/example.jpg)  

<h2 id="7">Tips for visualization</h2>

A plasmid visualization component using AngularJS is integrated into our pipeline, all the plasmid maps are in HTML formats. In order to view the map locally, you need to download the `js` folder in addtion to the HTML files.  

<h2 id="8">Citation</h2>  
Coming soon!  

<h2 id="8">Contact</h2>  
You Che, hkucheyou@gmail.com 
