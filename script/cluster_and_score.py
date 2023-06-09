import sys
import shutil
import os
import itertools
import multiprocessing
import subprocess
import numpy as np

import pycharmm
import pycharmm.lib as lib
import pycharmm.psf as psf
import pycharmm.read as read
import pycharmm.lingo as lingo
import pycharmm.energy as energy
import pycharmm.minimize as minimize
import pycharmm.generate as gen
import pycharmm.settings as settings
import pycharmm.grid as grid
import pycharmm.coor as coor
import pycharmm.write as write

# find the current working directory to use below
cwd = os.getcwd()
sys.path.append(f'{cwd}/script')
from FASTDock import dock_probe, cluster_n_score, find_rep, get_cluster, score

## Set FASTDock parameters
exec(open(f'{cwd}/script/parameters.py').read())

# Establish the directory for the test run
testdir = f'{cwd}/{testdir}'
if not os.path.exists(testdir): 
    print('No test directory {testdir}. Need to run script/dock_probes.py first')
    exit()
os.chdir(testdir)
if shutil.which('convpdb.pl') == None:
    print('Cannot find MMTSB Tool convpdb.pl, need MMTSB Tools installed')
    exit()

## Read in the topology and parameter file 
bl = settings.set_bomb_level(-1)
wl = settings.set_warn_level(-5)
v = settings.set_verbosity(-5)
read.rtf(f'"{cwd}/data/toppar/top_all36_prot.rtf"')
read.rtf(f'"{cwd}/data/toppar/top_all36_cgenff.rtf"', append = True)
read.rtf(f'"{cwd}/data/toppar/top_all36_cgenff.rtf"', append = True)
read.rtf(f'"{cwd}/data/toppar/probes.rtf"', append = True)
read.prm(f'"{cwd}/data/toppar/par_all36m_prot.prm"', flex = True)
read.prm(f'"{cwd}/data/toppar/par_all36_cgenff.prm"', append = True, flex = True)
read.prm(f'"{cwd}/data/toppar/probes.prm"', append = True, flex = True)
settings.set_bomb_level(bl)
settings.set_warn_level(wl)
settings.set_verbosity(v)

###############################
# #       CLUSTER PROBES      ##
# ##############################

# Cluster and score chemotypes
if os.path.exists('clusters'):
        shutil.rmtree('clusters')
        if os.path.exists('cluster-n-score.log'): os.remove('cluster-n-score.log')
        if os.path.exists('top.txt'): os.remove('top.txt')

os.makedirs('clusters/reps')

pool = multiprocessing.Pool()
pool.map(cluster_n_score, list(zip(itertools.repeat(cwd),fragnames)))
pool.close()

# get cluster representatives
for name in fragnames:
    j = fragnames.index(name)
    j = j*100
    find_rep(name, nsave, j)

# Check whether file exists so we don't append to inadvertent file
if os.path.exists('before_clustering.pdb'): os.remove('before_clustering.pdb')
# generate one file with final poses
file = open('top.txt', 'r')
files = []
for line in file.readlines():
    files.append(line.rstrip())
file.close()
for name in sorted(files):
    with open('before_clustering.pdb', 'ab') as outFile:
        with open(f'clusters/reps/{name}', 'rb') as com:
            outFile.write(com.read())

os.system('sed -i /REMARK/d before_clustering.pdb')
os.system('sed -i /END/d before_clustering.pdb')
os.system('convpdb.pl -segnames -renumber 1 before_clustering.pdb > tmp')
os.system('sed -i s/PRO0/LIG/g tmp')
os.system('mv tmp before_clustering.pdb')
os.system('rm top.txt')

#######################################
# #  GEN FASTDock CLUSTERS AND SCORE  ##
# ######################################

## List highest scoring pose for each probe type
top = []
for name in fragnames:
    # Sort all entries for {name} based on field 12 (cluster size)
    command = f"grep {name} cluster-n-score.log | sort -grk 12 "
    # Advance cluster only if cluster membership is >= {minsize}
    command += f"| awk ' $NF >= {minsize} ' | head -n 1 | awk '{{print $6}}'"
    pose = subprocess.check_output(command, shell=True)
    pose = pose.strip().decode('utf-8')
    top.append(f'{name}-{pose}')

if os.path.exists('final'):
        shutil.rmtree('final')
        if os.path.exists('score.txt'): os.remove('score.txt')

os.mkdir('final')

for probe in top:
    get_cluster(probe, top)

# get uniqe clusters and score
os.chdir('final')

os.mkdir('tmp')
for probe in top:
    # print cluster coordinates to tmp file
    os.system(f"grep ATOM {probe}_final.pdb | awk '{{print $6, $7, $8}}' > tmp/{probe}") 

os.chdir('tmp')
# check if coordinates for any clusters match
os.system("md5sum * | sort -u -k1,1 | awk '{print $2}' > scorelist") 
unique_clusters = np.loadtxt('scorelist',dtype=str)

os.chdir('../')
shutil.rmtree('tmp')
os.chdir('../')

for probe in unique_clusters: 
    score(probe,minsize,protein_psf,protein_pdb)
