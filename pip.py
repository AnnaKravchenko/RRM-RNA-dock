##!/usr/bin/env python3 

import os 
import sys
import argparse
import requests
import json
import numpy as np
import traceback
import random
import configparser

try:
    
    ATTRACTDIR = os.environ['ATTRACTDIR'] + '/'
except:
    print("Define the ATTRACTDIR environment variable with:")
    print("export ATTRACTDIR={your path to attract/bin}")
    print("or add that line in your /home/.bashrc")
    sys.exit()
try: 
    ATTRACTTOOLS = os.environ['ATTRACTTOOLS'] + '/'
except:
    print("Define the ATTRACTTOOLS environment variable with:")
    print("export ATTRACTTOOLS={your path to attract/tools}")
    print("or add that line in your /home/.bashrc")
    sys.exit()
try: 
    scr = os.environ['RRDOCK'] + '/'
except:
    print("Define the RRDOCK environment variable with:")
    print("export SRRDOCKCR={your path to RRM-RNA-dock/}")
    print("or add that line in your /home/.bashrc")
    sys.exit()
try: 
    nalib_source = os.environ['LIBRARY'] + '/'
except:
    print("Define the LIBRARY environment variable with:")
    print("export LIBRARY={your path to fraglib/}")
    print("or add that line in your /home/.bashrc")
    sys.exit()

def read_pdb_models(pdb):
    '''This function reads pdb file.
    If it contains more than 1 model it returns structure like [ [model 1] ... [model n] ]
    If it contains only one (or zero*) model it returns [ model 1 ].
    *Number of models is determined by the mumber of lines starting with "ENDMDL"  
    '''
    lines0 = open(pdb).readlines()
    lines, atoms, extralines, residues = [], [], [], []
    lines_up, atoms_up, residues_up = [], [], []
    flag = 0
    for l in lines0:
        if l.startswith("ENDMDL"):
            flag+=1
            
            atoms_up.append(np.array(atoms))
            lines_up.append(lines)
            residues_up.append(residues)
            lines, atoms, residues = [], [], []
        if not l.startswith("ATOM"):
            extralines.append((len(lines), l))
            continue
        x = float(l[30:38])
        y = float(l[38:46])
        z = float(l[46:54])
        i = int(l[23:26])

        atoms.append((x,y,z))
        lines.append(l)
        residues.append(i)
    if flag == 0:
        atoms_up = atoms 
        lines_up = lines 
        residues_up = residues
    if flag == 1:
        atoms_up = atoms_up[0]
        lines_up = lines_up[0] 
        residues_up = residues_up[0]

    return lines_up, np.array(atoms_up), extralines, np.array(residues_up)

def make_motif_boungfrag(rna_seq):
    motif = []
    bound = []
    frNum = len(rna_seq) - 2
    for n in range(frNum):
        head = n
        tail = n+3
        frag = rna_seq[head:tail]
        motif.append(frag)
        bound.append([n+1, frag])
    return motif, bound

def count_lines_for_restraints(fragIds, target, rnaSeq):
    flag = False
    x = 0
    if rnaSeq[fragIds[0]-1] in {'A','G'}:
        nuc1 = 7 
    else: 
        nuc1 = 6
    if rnaSeq[fragIds[1]-1] in {'A','G'}:
        nuc2 = 7
    else:
        nuc2 = 6

    if target == fragIds[0]:
            x = 0
    elif target == fragIds[1]:
            x = nuc1 
    elif target == fragIds[2]:
            x = nuc1 + nuc2
    
    if rnaSeq[target - 1 ] in {'A','G'}:
        flag = True 

    return x, flag 

def get_residue_id_for_first_ghost(workDir, proteinFile):
    fileName = "%s/%s" % (workDir, proteinFile) 
    with open(fileName, "r") as file:
        last_line = file.readlines()[-1]
    res = int(last_line[22:26])+1
    insert = "%4s" % res
    return insert

def get_reduced_nucleotides(filePathName):
    reducedNucleotides = {'RU', 'RA', 'RC', 'RG'}
    rnaLines = []
    lines, coord, extra, residues = read_pdb_models(filePathName)
    for i in range(len(lines)):
        if lines[i][18:20].strip() in reducedNucleotides:
            rnaLines.append(lines[i])
    return rnaLines

def rm_dir_precaution(path):
    flag = True
    print("Directory %s exists" % path)
    print("In order to proceed, it has to be empty")
    print("Empty this directory?")
    while flag:
        y = str(input('[y/n]: '))
        if y == 'n' or y =='y':
            flag = False 
    if y == 'n':
        sys.exit()
    else:
        cmd= "rm -r %s/* " % path 
        os.system(cmd)
    return 

def start_docking_precaution():
    flag = True
    print("Would you like to start the docking?")
    while flag:
        y = str(input('[y/n]: '))
        if y == 'n' or y =='y':
            flag = False 
    if y == 'n':
        flagOut = False 
    else:
        flagOut = True
    return flagOut

parser = argparse.ArgumentParser(description=__doc__)

parser.add_argument('-wdir', '--work_directory', type=str, required=True, 
help=""" Path to directory where RRM folder for the docking will be made.
""")

parser.add_argument('-id', '--uniProtID', type=str, required=True, default='P19339',
help=""" UniProt ID for the RRM, e.g. 'P31946' or 'P62258' etc.
""")

parser.add_argument('-rrm', '--rrm_domain_id', type=int, required=True, #nargs='+', 
help= """
RRM domain index e.g. '1' or '2' etc. 
""")

parser.add_argument('-seq', "--ss_rna_sequence", type=str, required=True, 
help= """
Single-stranded RNA sequence to be docked onto the RRM, at least 3 nucleotides long, e.g. 'CAC' or 'GCAC' etc.
""")
# TODO: make these 2 in one arg 2 position long. Or maybe not. Probably better not, it's ok like so.
parser.add_argument('-ancNucB1', "--anchoring_nucleotide_id_beta1",type=int, required=True, #nargs='+', 
help="""
Anchoring nucleotide index for betasheet 1, e.g. '1' or '2' etc. 
""")

parser.add_argument('-ancNucB3', '--anchoring_nucleotide_id_beta3', type=int, required=True, # nargs='+',  
help="""
Anchoring nucleotide index for betasheet 3, e.g. '2' or '3' etc.
""")

parser.add_argument('-conf', '--configuraion_file', type=str, const='/data3/akravche/data21/anchors/mdir/deliv1/config.ini', nargs='?',
help=""" Configuration file with advanced parameters. Leave empty to use default file. 
""")

#p.add_argument("--mutations", default=None, help="allow mutations on the protein")

args = parser.parse_args()

print("****************************************")
print("* Read user input...")
print("****************************************")
conf = str(args.configuraion_file)
if conf == "None":
    conf = '/data3/akravche/data21/anchors/mdir/deliv1/config.ini'

wdir = args.work_directory+'/'
prId = args.uniProtID
rrm = args.rrm_domain_id
rnaSeq =  args.ss_rna_sequence
nucB1 = args.anchoring_nucleotide_id_beta1
nucB3 = args.anchoring_nucleotide_id_beta3

# read from config file 
#config = ConfigParser()
#config = ConfigParser.ConfigParser()
#[]

try:
    config = configparser.ConfigParser()
    config.read(conf)
    mdir = config.get('DEFAULT', 'mdir')
    jsonFile = config.get("DEFAULT", "jsonFile")
    refB1 = config.get("DEFAULT", "refB1")
    refB3 = config.get("DEFAULT", "refB3")
    dist_restraint_B1 = config.get("DEFAULT", "dist_restraint_B1")
    dist_restraint_B3 = config.get("DEFAULT", "dist_restraint_B3")
    score_penalty = config.get("DEFAULT", "score_penalty")
except:
    print('Error while reading configuration file. ')
    print('Configuration file must have:')
    print(' - extention .ini  ')
    print(' - first line "[DEFAULT]"')
    print(' - all fields that are in original confiuraion file')
    sys.exit()

data = json.load(open(mdir+jsonFile))
proteins = []
for prot in data:
    proteins.append(prot)
    
if prId in proteins:
    #print("\033[1;32m This text is Bright Green \n")
    #print("\033[1;32m Valid UniProt index", prId)
    print("Valid UniProt index", prId)
else:
    #print("\033[1;31m *Invalid UniProt index: given -id/--uniProtID value (%s) is not in InteR3M database." % prId )
    #print("\033[1;31m Please re-run the script with valid UniProt index")
    print("*Invalid UniProt index: given -id/--uniProtID value (%s) is not in InteR3M database." % prId )
    print("Please re-run the script with valid UniProt index")
    sys.exit()

domains = []
for domain in data[prId]:
    domains.append(int(domain[-1]))

#for dom in args.rrm_domain_id:
if args.rrm_domain_id in domains: 
    print('Valid RRM domain index', rrm) 
else:
    print('*Invalid RRM domain index: given -rrm/--rrm_domain_id value (%s) is not in InteR3M database' % rrm)
    print('Please re-run the script with valid RRM domain index')
    print('Awailable domain indexes for %s : %s' % (prId, str(domains)[1:-1]) )
    sys.exit()

nucs='ACGU'
if len(rnaSeq) < 3:
    print("*Invalid ssRNA sequence: given -seq/--ss_rna_sequence (%s) contains less than 3 nucleotides" % rna)
    print('Please re-run the script with valid ssRNA sequence')
    sys.exit()
else: 
    if all(char in nucs for char in rnaSeq):
        print('Valid ssRNA sequence', rnaSeq)
    else: 
        print("*Invalid ssRNA sequence: given -seq/--ss_rna_sequence (%s) contains elements other then 'A', 'U', 'C', or 'G'" % rna)
        print('Please re-run the script with valid ssRNA sequence')
        sys.exit()

ancFlag = True 
if nucB1 > 0 :
    if nucB3 > 0:
        if nucB3 != nucB1:
            if len(rnaSeq) >= nucB3:
                if len(rnaSeq) >= nucB1:
                    print("Valid anchoring nucleotide index for Beta1", nucB1)
                    print("Valid anchoring nucleotide index for Beta3", nucB3)
                    spaceBtwnAnchors = abs(nucB1 - nucB3)
                    ancFlag = False 
                else:
                    print("*Invalid anchoring nucleotide index for Beta1: given -ancNucB1/--anchoring_nucleotide_id_beta1 value (%s) is out of range of given ssRNA sequence" % nucB1 )
            else:
                print("*Invalid anchoring nucleotide index for Beta3: given -ancNucB3/--anchoring_nucleotide_id_beta1 value (%s) is out of range of given ssRNA sequence" % nucB3 )
        else: 
            print("*Invalid anchoring nucleotides indexs: anchoring nucleotide for beta1 cannot be the same as  anchoring nucleotide for beta3")
    else:
        print("*Invalid anchoring nucleotide index for Beta3: index have to be positive integer")            
else:
    print("*Invalid anchoring nucleotide index for Beta1: index have to be positive integer") 
if ancFlag:
        print("Please re-run the script with valid anchoring nucleotide(s) index/indices")
        sys.exit()

fitting_region, domain = [] , [] 
fitting_region = data[prId]['RRM'+str(rrm)]['fitting_region']
domain  = data[prId]['RRM'+str(rrm)]['domain']

# print("********************************************************************************")
# print("* Load references ")
# print("********************************************************************************")

modelsB1 = read_pdb_models(mdir+refB1)[0] # now this is 5 models 
modelsB3 = read_pdb_models(mdir+refB3)[0] # now this is 5 models

numOfPrototypesB1= len(modelsB1) 
numOfPrototypesB3= len(modelsB3) 

print("********************************************************************************")
print("* Load protein 3D prediction from AlphaFold ")
print("********************************************************************************")

pdb="%s.pdb" % prId
if not os.path.isdir(wdir+prId):
    os.mkdir(wdir+prId)
rrmDir = wdir+prId+'/rrm'+str(rrm)
if not os.path.isdir(rrmDir):
    os.mkdir(rrmDir)
else:
    if len(os.listdir(rrmDir)) != 0 :
        # rm inside of rrmDir or exit script
        rm_dir_precaution(rrmDir)

#if not os.path.exists(rrmDir+'/proteinAFold.pdb'): 
URL = 'https://alphafold.ebi.ac.uk/files/AF-%s-F1-model_v2.pdb' % prId 
response = requests.get(URL)
open(rrmDir+'/proteinAFold.pdb', "wb").write(response.content)

print("********************************************************************************")
print("* Convert protein structure to sequence ")
print("********************************************************************************")

cmd = "python3 %s/anchoring/pdb2fasta.py %s > %s " % (scr, rrmDir+'/proteinAFold.pdb', rrmDir+'/proteinAFold.fasta')
os.system(cmd)

print("********************************************************************************")
print("* Extract & save domain of interest from protein structure ")
print("********************************************************************************")

lines, coord, extra, residues = read_pdb_models(rrmDir+'/proteinAFold.pdb')
domain_lines , domain_coord, = [],[] 

for i in range(domain[0], domain[-1]+1):
    for x in range(len(residues)):
        if residues[x] == i:
            domain_lines.append(lines[x])
with open(rrmDir+'/extract_domain.pdb', "w") as file:
    for el in domain_lines:
        file.write(el)

print("********************************************************************************")
print("* Renumber atoms in the domain ")
print("********************************************************************************")

cmd = "python3 %s/anchoring/renumber_atoms.py %s %s " % (scr,rrmDir+'/extract_domain.pdb', rrmDir+'/domain.pdb')
os.system(cmd)

lines, coord, extra, residues = [],[],[],[]
lines, coord, extra, residues = read_pdb_models(rrmDir+'/domain.pdb')

print("********************************************************************************")
print("* Get first fitting region from the domain ")
print("********************************************************************************")

fit1_lines = []
for i in range(fitting_region[0][0], fitting_region[0][-1]+1):
    for x in range(len(residues)):
        if residues[x] == i:
            fit1_lines.append(lines[x])

print("********************************************************************************")
print("* Get second fitting region from the domain ")
print("********************************************************************************")

fit3_lines = []
for i in range(fitting_region[1][0], fitting_region[1][-1]+1):
    for x in range(len(residues)):
        if residues[x] == i:
            fit3_lines.append(lines[x])
            
print("********************************************************************************")
print("* Softlink the references and split each in the models  ")
print("********************************************************************************")
src1 = mdir+refB1
src2 = mdir+refB3
dst1 = rrmDir+'/protoB1.pdb'
dst2 = rrmDir+'/protoB3.pdb'

if not os.path.isfile(dst1):
    os.symlink(src1, dst1)
if not os.path.isfile(dst2):
    os.symlink(src2, dst2)

cmd= "python2 %s/splitmodel.py %s > /tmp/split.out " % (ATTRACTTOOLS, dst1 )
os.system(cmd)
cmd= "python2 %s/splitmodel.py %s > /tmp/split.out " % (ATTRACTTOOLS, dst2 )
os.system(cmd)

print("********************************************************************************")
print("* Renumber first model prototypes B1 & B3, to get reference atom numbers for fitting  ")
print("********************************************************************************")
i=1
protoName = rrmDir+'/protoB1-%s.pdb' % i
protoName_remunber =  rrmDir+'/protoB1-%s-renumber.pdb' % i
cmd = "python3 %s/anchoring/renumber_atoms.py %s %s " % (scr, protoName, protoName_remunber)
os.system(cmd)

linesRefB1, coordRefB1, extraRefB1, residuesRefB1 = read_pdb_models(protoName_remunber)

protoName = rrmDir+'/protoB3-%s.pdb' % i
protoName_remunber =  rrmDir+'/protoB3-%s-renumber.pdb' % i
cmd = "python3 %s/anchoring/renumber_atoms.py %s %s " % (scr, protoName, protoName_remunber)
os.system(cmd)

linesRefB3, coordRefB3, extraRefB3, residuesRefB3 = read_pdb_models(protoName_remunber)

print("********************************************************************************")
print("* Select atoms for the reference fitting in the first fitting region  ")
print("********************************************************************************")
# line[13:16].strip() is atom names, like in set2fit 
set2fit={'N','CA','C','O'} # backbone only
domainAtoms2fitB1 = []
for line in fit1_lines:
    if line[13:16].strip() in set2fit:
        domainAtoms2fitB1.append(int(line[7:11]))

refB1Atoms2fit = []
for line in linesRefB1:
    if line[13:16].strip() in set2fit:
        refB1Atoms2fit.append(int(line[7:11]))

print("********************************************************************************")
print("* Select atoms for the reference fitting in the second fitting region  ")
print("********************************************************************************")
domainAtoms2fitB3 = []
for line in fit3_lines:
    if line[13:16].strip() in set2fit:
        domainAtoms2fitB3.append(int(line[7:11]))

refB3Atoms2fit = []
for line in linesRefB3:
    if line[13:16].strip() in set2fit:
        refB3Atoms2fit.append(int(line[7:11]))

print("********************************************************************************")
print("* Fit all beta1 reference models to selections from the fitting regions ")
print("********************************************************************************")

for i in range(1,numOfPrototypesB1+1):
    protoName = rrmDir+'/protoB1-%s.pdb' % i
    outName = rrmDir+'/fitted-protoB1-%s.pdb' % i
    cmd= "python3 %s/fit.py %s %s --selection1 %s --selection2 %s > %s " % (ATTRACTTOOLS, rrmDir+'/domain.pdb' , protoName , str(domainAtoms2fitB1)[1:-1].replace(',',''), str(refB1Atoms2fit)[1:-1].replace(',',''), outName )
    os.system(cmd)

for i in range(1,numOfPrototypesB3+1):
    protoName = rrmDir+'/protoB3-%s.pdb' % i
    outName = rrmDir+'/fitted-protoB3-%s.pdb' % i
    cmd= "python3 %s/fit.py %s %s --selection1 %s --selection2 %s > %s " % (ATTRACTTOOLS, rrmDir+'/domain.pdb' , protoName , str(domainAtoms2fitB3)[1:-1].replace(',',''), str(refB3Atoms2fit)[1:-1].replace(',',''), outName )
    os.system(cmd)

print("********************************************************************************")
print("* Reduce domain ")
print("********************************************************************************")
cmd = 'python  %s/reduce.py %s >> /tmp/reduce.out ' % (ATTRACTTOOLS , rrmDir+'/domain.pdb')
os.system(cmd)

lenProtr = sum(1 for line in open(rrmDir+'/domainr.pdb'))
#print(lenProtr)
# HERE ITS GONNA BE DIFFERENT REGARDING OF THE NUMBER OF FRAGMENTS TO DOCK 

if spaceBtwnAnchors <= 2:
    #print('ancNucs are in one fragment')
    #same from another 
    print("********************************************************************************")
    print("* Get just nucleotides from fitted prototypes ")
    print("********************************************************************************")
    
    rnasB1=[]
    for i in range(1, numOfPrototypesB1+1):
        protoName = rrmDir+'/fitted-protoB1-%s.pdb' % i
        reducedRNA = get_reduced_nucleotides(protoName)
        rnasB1.append(reducedRNA)

    rnasB3=[]
    for i in range(1, numOfPrototypesB3+1):
        protoName = rrmDir+'/fitted-protoB3-%s.pdb' % i
        reducedRNA = get_reduced_nucleotides(protoName)
        rnasB3.append(reducedRNA) 
    # not same 
    print("********************************************************************************")
    print("* Create dir structure, create reseptor file with 2 ghost nucleotides  ")
    print("********************************************************************************")
    ins = get_residue_id_for_first_ghost(rrmDir, 'domainr.pdb')
    domainOrig = read_pdb_models(rrmDir+'/domainr.pdb')[0]
    for i in range(1, numOfPrototypesB1+1):
        domain2use = domainOrig[:]
        nucleotide = rnasB1[i-1][:]
        for l in nucleotide:
            domain2use.append(l[:22]+ins+l[26:56]+' 99 '+l[60:])
        pathDir1 = rrmDir+'/b1_%s' % i 
        if not os.path.isdir(pathDir1):
            os.mkdir(pathDir1)
        # here I should have domainr+1ghost
        for j in range(1, numOfPrototypesB3+1):
            domain2use_secongGost = domain2use[:]
            ins2 = "%4s" % str(int(ins)+1)
            nucleotide = rnasB3[j-1][:]
            for l in nucleotide:
                domain2use_secongGost.append(l[:22]+ins2+l[26:56]+' 99 '+l[60:])
            pathDir2 = pathDir1+'/b3_%s' % j
            if not os.path.isdir(pathDir2):
                os.mkdir(pathDir2)

            dst = pathDir2+'/domainr-protoB1_%s-B3_%s-tmp.pdb' % (i , j)
            with open(dst, "w") as file:
                for el in domain2use_secongGost:
                    file.write(el)
            #domain2use = []
            if not os.path.isdir(pathDir2+'/nalib'):
                os.symlink(nalib_source, pathDir2+'/nalib')
            dst2 = pathDir2+'/domainr-protoB1_%s-B3_%s.pdb' % (i , j)
            cmd = "python3 %s/anchoring/renumber_atoms.py %s %s" % (scr, dst, dst2)
            os.system(cmd)
            os.remove(dst)

    print("********************************************************************************")
    print("* Work out which fragment will be docked  ")
    print("********************************************************************************")
    
    motif, bound = make_motif_boungfrag(rnaSeq)
    with open(rrmDir+'/motif.list', "w") as file:
        for el in list(set(motif)):
            file.write(el+'\n')

    with open(rrmDir+'/boundfrag.list', "w") as file:
        for el in bound:
            file.write("%s %s\n" % (el[0], el[1]))
    frag = []
    if nucB1 < nucB3: 
        n1 = nucB1
        n2 = nucB3
    else:
        n1 = nucB3
        n2 = nucB1
    if spaceBtwnAnchors == 1: 
        leftSide = len(rnaSeq[:n1-1])
        rightSide = len(rnaSeq[n2:])
        if leftSide < rightSide:
            frag = [n1,n2,n2+1]
        elif leftSide > rightSide:
            frag = [n1-1,n1,n2]
        elif leftSide == rightSide:
            n3 =  random.sample([n1-1, n2+1], 1)[0]
            if n3 == n1-1:
                frag = [n1-1, n1,n2]
            else:
                frag = [n1,n2,n2+1]
    elif spaceBtwnAnchors == 2:
        frag = [n1,n1+1,n2]
    
    fragSequence =  rnaSeq[frag[0]-1:frag[2]]
    with open(rrmDir+'/frag.info', "w") as file:
        for el in [frag, nucB1, nucB3, rnaSeq[frag[0]-1:frag[2]]]:
            file.write(str(el)+'\n')

    print("********************************************************************************")
    print("* Create restraint files  ")
    print("********************************************************************************")
    # i have really fuzzy understanding of how it works, it really not a good code.
    # but I'm one weekend before writting time and it's working so just let it be 

    moveGhost = 0
    ghostOrder = 1
    restrFile_p1, restrFile_p2, finRest = [], [], []
    for nucleotide in [n1,n2]:
        if nucleotide == nucB1:
            anc = 'B1'
            lenGhostRNA = len(rnasB1[0])
            dist_restraint = dist_restraint_B1
        else: 
            anc = 'B3'
            lenGhostRNA = len(rnasB3[0])
            dist_restraint = dist_restraint_B3
        j = 0
        restrFile = []
        x, skipFirstNucOfBase = count_lines_for_restraints(frag, nucleotide, rnaSeq)

        for i in range(1,lenGhostRNA+1):
            if skipFirstNucOfBase:
                if i >= 4:
                    j=1
            ghostLine = 'ghost_B%s-%s 1 %s' % (ghostOrder, i,lenProtr+i+moveGhost)
            realLine = 'real_%s-%s 1 %s' % (anc, i,lenProtr+i+lenGhostRNA*2+x+j)        
            restrFile.append(ghostLine)
            restrFile.append(realLine)
        restrFile_p1 = restrFile_p1 + (restrFile)
        restrFile = []
        
        for i in range(1,lenGhostRNA+1):
            restrFile_p2.append("ghost_%s-%s real_%s-%s 1 %s %s" % (anc, i, anc,i, dist_restraint, score_penalty))
        moveGhost+= 6 
        ghostOrder = 3

    finRest = restrFile_p1[:]
    finRest = finRest + [' ']
    finRest = finRest + restrFile_p2

    for b1 in range(1,len(rnasB1)+1):
        for b3 in  range(1,len(rnasB3)+1):
            p = (rrmDir+'/b1_%s/b3_%s/' % (b1,b3) )
            with open(p+'/restraints.txt', "w") as file:
                for el in finRest:
                    file.write(el+'\n')

    with open(rrmDir+'/restraints.txt', "w") as file:
        for el in finRest:
            file.write(el+'\n')

else:

    print("********************************************************************************")
    print("* Get just nucleotides from fitted prototypes ")
    print("********************************************************************************")
    
    def get_reduced_nucleotides(filePathName):
        reducedNucleotides = {'RU', 'RA', 'RC', 'RG'}
        rnaLines = []
        lines, coord, extra, residues = read_pdb_models(filePathName)
        for i in range(len(lines)):
            if lines[i][18:20].strip() in reducedNucleotides:
                rnaLines.append(lines[i])
        return rnaLines

    rnasB1=[]
    for i in range(1, numOfPrototypesB1+1):
        protoName = rrmDir+'/fitted-protoB1-%s.pdb' % i
        reducedRNA = get_reduced_nucleotides(protoName)
        rnasB1.append(reducedRNA)

    rnasB3=[]
    for i in range(1, numOfPrototypesB3+1):
        protoName = rrmDir+'/fitted-protoB3-%s.pdb' % i
        reducedRNA = get_reduced_nucleotides(protoName)
        rnasB3.append(reducedRNA)   
    
    print("********************************************************************************")
    print("* Change nucleotide residue number, atom types (to ghost) and   ")
    print("* add it to the end of reduced domain file  ")
    print("********************************************************************************")  
       
    ins = get_residue_id_for_first_ghost(rrmDir, 'domainr.pdb')
    
    domainr_body_orig = read_pdb_models(rrmDir+'/domainr.pdb')[0]

    for i in range(1,numOfPrototypesB1+1):
        domainr_body = domainr_body_orig[:]
        dst = rrmDir+'/domainr-protoB1-%s-tmp.pdb' % i
        rna = rnasB1[i-1][:]
        for l in rna:
            domainr_body.append(l[:22]+ins+l[26:56]+' 99 '+l[60:])
        with open(dst, "w") as file:
            for el in domainr_body:
                file.write(el)

    for i in range(1,numOfPrototypesB3+1):
        domainr_body = domainr_body_orig[:]
        dst = rrmDir+'/domainr-protoB3-%s-tmp.pdb' % i
        rna = rnasB3[i-1][:]
        for l in rna:
            domainr_body.append(l[:22]+ins+l[26:56]+' 99 '+l[60:])
        with open(dst, "w") as file:
            for el in domainr_body:
                file.write(el)

    print("********************************************************************************")
    print("* Renumber ghost atoms in the domain ")
    print("********************************************************************************")
    for i in range(1,numOfPrototypesB1+1):
        src = rrmDir+'/domainr-protoB1-%s-tmp.pdb' % i
        dst = rrmDir+'/domainr-protoB1-%s-done.pdb' % i
        cmd = "python3 %s/anchoring/renumber_atoms.py %s %s" % (scr, src, dst)
        os.system(cmd)
    for i in range(1,numOfPrototypesB3+1):
        src = rrmDir+'/domainr-protoB3-%s-tmp.pdb' % i
        dst = rrmDir+'/domainr-protoB3-%s-done.pdb' % i
        cmd = "python3 %s/anchoring/renumber_atoms.py %s %s" % (scr, src, dst)
        os.system(cmd)

    print("********************************************************************************")
    print("* Remove used files  ")
    print("********************************************************************************")
    os.remove(rrmDir+'/extract_domain.pdb')
    os.remove(rrmDir+'/protoB1-1-renumber.pdb')
    os.remove(rrmDir+'/protoB3-1-renumber.pdb')

    for i in range(1,numOfPrototypesB1+1):
        filename1 = rrmDir+'/domainr-protoB1-%s-tmp.pdb' % i
        filename2 = rrmDir+'/protoB1-%s.pdb' % i 
        filename3 = rrmDir + '/fitted-protoB1-%s.pdb' % i
        os.remove(filename1)
        os.remove(filename2)
        os.remove(filename3)

    for i in range(1,numOfPrototypesB3+1):
        filename1 = rrmDir+'/domainr-protoB3-%s-tmp.pdb' % i
        filename2 = rrmDir+'/protoB3-%s.pdb' % i 
        filename3 = rrmDir + '/fitted-protoB3-%s.pdb' % i
        os.remove(filename1)
        os.remove(filename2)
        os.remove(filename3)  

    print("********************************************************************************")
    print("* Make dirs & move files to corresponding ones   ")
    print("********************************************************************************")
    dirsLevel1 = ['b1', 'b3']
    dirsLevel2 = [range(1,numOfPrototypesB1+1) , range(1,numOfPrototypesB3+1)]
    
    for i in range(len(dirsLevel1)):
        pathDir = rrmDir+'/%s' % dirsLevel1[i]
        x =dirsLevel1[i][1]
        if not os.path.isdir(pathDir):
            os.mkdir(pathDir)
        for j in range(len(dirsLevel2[i])):
            pathDir2 = pathDir+'/%s' % dirsLevel2[i][j]
            y=dirsLevel2[i][j]
            if not os.path.isdir(pathDir2):
                os.mkdir(pathDir2)
            if not os.path.isdir(pathDir2+'/nalib'):
                os.symlink(nalib_source, pathDir2+'/nalib')
            mvSrc = rrmDir+"/domainr-protoB%s-%s-done.pdb" % (x,y) 
            mvDst = pathDir2+"/domainr-protoB%s-%s-done.pdb" % (x,y)
            cmd = "mv %s %s" % (mvSrc, mvDst)
            os.system(cmd)

    print("********************************************************************************")
    print("* Work out which fragments will be docked  ")
    print("********************************************************************************")
    
    motif, bound = make_motif_boungfrag(rnaSeq)
    with open(rrmDir+'/motif.list', "w") as file:
        for el in list(set(motif)):
            file.write(el+'\n')

    with open(rrmDir+'/boundgrag.list', "w") as file:
        for el in bound:
            file.write("%s %s\n" % (el[0], el[1]))
    frag1, frag2 = [], []
    if nucB1 < nucB3: 
        n1 = nucB1
        n2 = nucB3
    else:
        n1 = nucB3
        n2 = nucB1

    if  (spaceBtwnAnchors >= 3 ) and (spaceBtwnAnchors <= 5):
        frag1 = [n1,n1+1,n1+2]
        frag2 = [n2-2,n2-1,n2]
    else:
        if n1-1 != 0:
            frag1 = [n1-1,n1,n1+1]
        else:
            frag1 = [n1,n1+1,n1+2]
        if n2 != len(rnaSeq):
            frag2 = [n2-1,n2,n2+1]
        else:
            frag2 = [n2-2,n2-1,n2]      
    print(frag1)
    print(frag2)

    print("********************************************************************************")
    print("* Create restraint files  ")
    print("********************************************************************************")

    #if spaceBtwnAnchors > 2: # then 2 fragments are docked separately 
    for nucleotide in [n1,n2]:
        if nucleotide in frag1:
            frag = frag1[:]
        else:
            frag = frag2[:]
        if nucleotide == nucB1:
            anc = 'B1'
            lenGhostRNA = len(rnasB1[0])
            dist_restraint = dist_restraint_B1
        else: 
            anc = 'B3'
            lenGhostRNA = len(rnasB3[0])
            dist_restraint = dist_restraint_B3
        j = 0
        restrFile = []
        x, skipFirstNucOfBase = count_lines_for_restraints(frag, nucleotide, rnaSeq)

        with open(rrmDir+'/frag%s.info' % anc, "w") as file:
            for el in [frag, nucleotide, rnaSeq[frag[0]-1:frag[2]]]:
                file.write(str(el)+'\n')

        for i in range(1,lenGhostRNA+1):
            if skipFirstNucOfBase:
                if i >= 4:
                    j=1
            ghostLine = 'ghost_%s 1 %s' % (i,lenProtr+i)
            realLine = 'real_%s 1 %s' % (i,lenProtr+i+lenGhostRNA+x+j)        
            restrFile.append(ghostLine)
            restrFile.append(realLine)
        restrFile.append('\n')
        for i in range(1,lenGhostRNA+1):
            restrFile.append("ghost_%s real_%s 1 %s %s" % (i, i, dist_restraint, score_penalty))

        fileName = rrmDir+'/restraint-%s.txt' % anc
        with open(fileName, "w") as file:
            for el in restrFile:
                file.write(el+'\n')
                
    print("********************************************************************************")
    print("* Move restraint files to corresponding directories ")
    print("********************************************************************************")

    for i in range(len(dirsLevel1)):
        pathDir = rrmDir+'/%s' % dirsLevel1[i]
        x =dirsLevel1[i][1]
        for j in range(len(dirsLevel2[i])):
            pathDir2 = pathDir+'/%s' % dirsLevel2[i][j]
            cpSrc = rrmDir+"/restraint-B%s.txt" % x 
            cpDst = pathDir2+"/restraint-B%s.txt" % x
            cmd = "cp %s %s" % (cpSrc, cpDst)
            os.system(cmd)


def dock(dockConfig,dockDir,motif,proteinr,restraints):
    try:
        RANDSEARCH = os.environ['RANDSEARCH'] + '/'
    except:
        print("Define the RANDSEARCH environment variable with:")
        print("export RANDSEARCH={your path to randsearch/}")
        print("or add that line in your /home/.bashrc")
        sys.exit()
    try:
        config = configparser.ConfigParser()
        config.read(dockConfig)
        cpu = config.get('DOCKING', 'docking_cpu')
        docking_tmp = config.get('DOCKING', 'docking_tmp')
        if docking_tmp == '0':
            docking_tmp = dockDir

    except:
        print('Error while reading configuration file. ')
        sys.exit()

    cmd = "bash %s/anchoring/dock_restr.sh %s %s %s %s %s %s"  % (scr, dockDir, motif, cpu, docking_tmp, proteinr, restraints) 
    #print("bash %s/anchoring/dock_restr.sh %s %s %s %s %s %s"  % (scr, dockDir, motif, cpu ,docking_tmp, proteinr, restraints) )
    os.system(cmd)
    print("********************************************************************************")
    print("* Docking for %s is done " % proteinr )
    print("********************************************************************************")

if start_docking_precaution() :

    for i in range(1,len(rnasB1)+1):
        for j in range(1,len(rnasB3)+1):
            
            wdirDock = "%s/b1_%s/b3_%s/" % (rrmDir,i,j)
            proteinr = 'domainr-protoB1_%s-B3_%s.pdb' % (i,j)  
            dock(conf,wdirDock,fragSequence,proteinr,'restraints.txt')