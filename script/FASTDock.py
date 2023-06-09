# Set of functions for use in the FASTDock workflow
# Written by Furyal Ahmed (1/23) 

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
import pycharmm.nbonds as nbonds
import pycharmm.charmm_file as charmm_file

import numpy as np
import pickle
import pandas as pd
import sys
import argparse
import random
import re
import subprocess
import os
import csv

orig_out = sys.stdout

def process_scores(name):

    col_names = ['nr','gridE','totE','gvdW','gElec']
    s = open(f'scores/dock{name}.dat', 'r')
    s_dict = {}
    for i in col_names:
        s_dict[i] = []
    for l in s:
        if len(l) == 0: continue
        sl = l.split()
        for i in range(len(col_names)):
            s_dict[col_names[i]].append(sl[i])
    s.close()
    return s_dict

def process_clusters(name):
    c = open(f'clusters/{name}/{name}.out', 'r')
    c_dict = {}
    r_dict = False
    for l in c:
        if l[0] == '#': continue
        elif l[0] == '@' and not l.split()[1] == 't':
            key = l.split()[1]
            c_dict[key] = []
            r_dict = True
        elif r_dict:
            c_dict[key].append(l.split())
        else: continue
    return c_dict

def cluster_n_score(args):
    cwd = args[0]
    i = args[1]
    kbt = 0.6 # kcal/mol
    os.system(f'mkdir clusters/{i}')
    os.system(f'cp {cwd}/data/probe_pdb/{i}.pdb clusters/{i}')
    wd = os.getcwd()
    os.chdir('clusters')

    cluster = 'cluster.pl -kclust -radius 4 -iterate -mode rmsd'
    cluster += ' -pdb -selmode heavy -nolsqfit -centroid '
    cluster += f'-centout {i}/{i} -log {i}/{i}.log'
    cluster += f' ../pdb_dock/dock_{i}_*.pdb > {i}/{i}.out'
    os.system(cluster)
    os.chdir(wd)
    score = process_scores(i) # returns a dictionary
    clusters = process_clusters(i) # returns a dictionary
    for k in clusters.keys():  
        q = np.exp([-kbt*float(score['gridE'][int(clusters[k][i][0])-1]) \
                    for i in range(len(clusters[k]))])
        e = [float(score['gridE'][int(clusters[k][i][0])-1]) \
             for i in range(len(clusters[k]))]
        eave = np.sum(np.multiply(e,q)) / np.sum(q)
        G_score =  eave -\
        kbt*float(len(clusters[k]))/float(len(score['gridE']))*\
        np.log(float(len(clusters[k]))/float(len(score['gridE'])))
        with open('cluster-n-score.log', 'a') as f:
                sys.stdout = f
                print(f'Average Score for cluster {i} {k} is {G_score} '\
                f'number of members {len(clusters[k])}')
                sys.stdout = orig_out
        with open('top.txt', 'a') as f:
                sys.stdout = f
                if len(clusters[k]) >= 10:
                        print(f'{i}-{k}.pdb')
                sys.stdout = orig_out


def find_rep(name,nsave, j):
        skip = nsave + 4 # skip mmtsb comments in file
        os.chdir(f'clusters/{name}')
        df = pd.read_csv(f'{name}.out',comment='@', dtype=str, header=None, skiprows=skip, sep=' ')
        data = df[[0,2]].to_numpy()

        os.system("grep '@' "+name+".out | tail -n +2 | awk '{print $4}' > count")
        frame_count = np.loadtxt("count",dtype=int)
        clusters = subprocess.check_output("ls *-t* | wc -l ", shell=True)
        clusters = int(clusters)
        temp = [0]*clusters # number of clusters in total
        
        if len(temp)==1:
                t1 = data[:,1]
                temp = min(t1)

                reps = []
                rep = df.loc[df[2] == f'{temp}'].to_numpy()[:,1]
                if rep.size != 0:
                      reps.append(rep[0])

        else:
                t1 = data[0:frame_count[0],1]
                temp[0] = min(t1)
        
                for i in range(0,clusters-1):
                    t2 = data[np.sum(frame_count[0:i+1]):np.sum(frame_count[0:i+2]),1]
                    a = min(t2)
                    temp[i+1] = a
       
                reps = []
                for i in range(0,clusters):
                        rep = df.loc[df[2] == f'{temp[i]}'].to_numpy()[:,1]
                        if rep.size != 0:
                              reps.append(rep[0])

        for i in range(0,len(reps)):
                os.system(f'convpdb.pl -segnames -renumber {j+i+1} ../{reps[i]} > tmp')
                os.system('sed -i s/PRO0/LIG/g tmp')
                os.system(f'mv tmp ../reps/{name}-t.{i+1}.pdb')

        os.chdir('../../')



def get_cluster(probe, top):
    os.system(f'cp clusters/reps/{probe}.pdb test.pdb')
    os.system('convpdb.pl -segnames -renumber 999 test.pdb > tmp')
    os.system('sed -i s/PRO0/LIGA/g tmp')
    os.system('mv tmp test.pdb')

    # read in probe for cluster center
    read.sequence_pdb('test.pdb')
    gen.new_segment(seg_name = 'LIGA')
    read.pdb('test.pdb', resid = True)

    # read in all docked probes
    read.sequence_pdb('before_clustering.pdb')
    gen.new_segment(seg_name = 'LIG')
    read.pdb('before_clustering.pdb', resid = True)

    ## create initial cluster based on
    ## all-atom distance in a 4A radius
    center = pycharmm.SelectAtoms(seg_id = 'LIGA')
    init_cluster = pycharmm.SelectAtoms.around(center, 4).whole_residues()#.by_seg_id("LIGA")

    # calculate geometric center of initial cluster
    stats = coor.stat(init_cluster)
    xcen = stats['xave'] 
    ycen = stats['yave']
    zcen = stats['zave']

    ## select probes 9A from center of 
    ## initial cluster to make final cluster
    index = top.index(probe)
    cluster = init_cluster.in_sphere(xcen, ycen, zcen, radius=9).whole_residues()\
            & ~pycharmm.SelectAtoms(seg_id='LIGA')
    write.coor_pdb(f'final/{probe}_final.pdb',selection=cluster)
    psf.delete_atoms(pycharmm.SelectAtoms().all_atoms())
    os.system('rm test.pdb')


## Rank probe clusters according to 
## ratio of contacts w/ protein to probes
def score(probe,minsize,protein_psf,protein_pdb):
    # read in protein
    read.psf_card(f'{protein_psf}', append = True)
    read.pdb(f'{protein_pdb}', resid = True)
    coor.orient(by_rms=True,by_mass=False,by_noro=False)
    
    no_probes = int(subprocess.check_output(f"grep ATOM final/{probe}_final.pdb | awk '{{print $5}}' | uniq | wc -l", shell=True))
    print(f'Found {no_probes} in {probe}')
    if no_probes > 0:
        # read in FASTDock cluster
        read.sequence_pdb(f'final/{probe}_final.pdb')
        gen.new_segment(seg_name = 'LIG')
        read.pdb(f'final/{probe}_final.pdb', resid = True)

        # calculate protein - ligand contacts
        lig = pycharmm.SelectAtoms(seg_id='LIG')
        contacts = pycharmm.SelectAtoms.around(lig, 4) \
                & ~pycharmm.SelectAtoms(seg_id='LIG')
        write.coor_pdb('tmp.pdb',selection=contacts)

        no_probes = int(subprocess.check_output(f"grep ATOM final/{probe}_final.pdb | awk '{{print $5}}' | uniq | wc -l", shell=True))
        no_contacts = int(subprocess.check_output(f"grep ATOM tmp.pdb | awk '{{print $5}}' | wc -l",shell=True))
        ratio = no_contacts / no_probes

        if no_probes >= minsize:
            with open('score.txt', 'a') as f:
                sys.stdout = f
                print('{:0.3f}'.format(ratio), probe)
                sys.stdout = orig_out
    else:
        print(f'Found probes for {probe}, skipping')
    psf.delete_atoms(pycharmm.SelectAtoms().all_atoms())
    os.system('rm tmp.pdb')


## This script is to dock multiple probes to
## protein grid using FFTDock
def dock_probe(cwd='./',addmin = True, randlig = True, fragnames = [], nsave = 2000, unit = 30):
    '''
    fragnames: array with names of probe pdbs
    nsave: number of poses to save per probe
    unit: unit number that the protein grid is read into
    '''
    for name in fragnames:
        read.psf_card(f'"{cwd}/data/probe_pdb/{name}.psf"') #embed in double-quotes to
        read.pdb(f'"{cwd}/data/probe_pdb/{name}.pdb"')      # ensure case retained

        # read the grid onto the GPU
        grid = pycharmm.CharmmFile(file_name='grid/protein_grid.bin',\
                                file_unit=unit,formatted=False, read_only=True)
        grid.open()
        lingo.charmm_script(f'fftg read unit {grid.file_unit}')
        grid.close()
        
        # Read in quaternions to sample rotation space of fragment
        qua = pycharmm.CharmmFile(file_name=f'"{cwd}/data/toppar/fftdock_rotation_2.qua"',\
                                file_unit=unit+1,formatted=True, read_only=True)
        qua.open()
        if randlig:
            lingo.charmm_script(f'''
            fftg lcon ncon 1 icon 1 nrok {nsave} -
            sizb 100 rnqu quau {qua.file_unit} select segid LIG end
            ''')
        else:
            lingo.charmm_script(f'''
            fftg lcon ncon 1 icon 1 nrok {nsave} -
            sizb 100 quau {qua.file_unit} select segid LIG end
            ''')
        qua.close()

        grid.open()
        lingo.charmm_script(f'grid read unit {grid.file_unit} select segid lig end')
        grid.close()
       
        scores = []        
        for nr in range(1,nsave+1):
            lingo.charmm_script(f'fftg coor icon 1 irot {nr}')
            if addmin: minimize.run_sd(nstep=100)
            energy.show()
            ener = energy.get_total()
            grvdW = lingo.get_energy_value('GRVD')
            grelec = lingo.get_energy_value('GREL')
            egrid = grvdW + grelec
            write.coor_pdb(f'pdb_dock/dock_{name}_{nr}.pdb',
                           title=f'* Solution {nr} {egrid} {ener} {grvdW} {grelec}')
            scores.append([nr, egrid, ener, grvdW, grelec])

        with open(f'scores/dock{name}.dat', "w") as f:
            csv.writer(f, delimiter=' ').writerows(scores)

        lingo.charmm_script("grid clear")
        lingo.charmm_script("fftg clear")
        psf.delete_atoms(pycharmm.SelectAtoms().all_atoms())         
    return
