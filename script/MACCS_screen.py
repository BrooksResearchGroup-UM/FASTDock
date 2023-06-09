from rdkit import Chem
from rdkit.Chem import ChemicalFeatures
from rdkit import RDConfig
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.Chem.Pharm2D.SigFactory import SigFactory
from rdkit.Chem.Pharm2D import Generate
from rdkit.Chem.Pharm2D import Gobbi_Pharm2D
from rdkit.Chem import MACCSkeys

import os
import shutil
import csv
import sys 
import numpy as np

# find the current working directory to use below
cwd = os.getcwd()
sys.path.append(f'{cwd}/script')

## Set FASTDock parameters
exec(open(f'{cwd}/script/parameters.py').read())

# Establish the directory for the test run
testdir = f'{cwd}/{testdir}'
if not os.path.exists(testdir): 
    print(
        f'''Missing test directory {testdir}, 
        have you run dock_probs and cluster_and_score?
        If not, run these first''')
    exit()
os.chdir(testdir)

########################
# #     PARAMETERS     ##
# #######################

# Read the scores from cluster_and_score.py output then pick cluster with largest score
s = open(f'{testdir}/score.txt','r')
scores = {}
for l in s:
    scores[l.split()[1]] = float(l.split()[0])
scores = dict(sorted(scores.items(), key=lambda item: item[1],reverse=True))
ref_cluster = list(scores.keys())[0]
refpdb = f'{testdir}/final/{ref_cluster}_final.pdb' # path to ligand file for reference molecule
database = f'{cwd}/data/database.sdf' # path to sdf file containing database for screening
filename = f'{testdir}/final/{ref_cluster}_tanimoto' # name for score output
ref = f'{testdir}/final/{ref_cluster}_final.mol2'
if shutil.which('obabel') == None:
    print('Cannot find openbabel obabel, need obabel installed')
    exit()
# Use openbabel to make refpdb into ref (mol2)
os.system(f'obabel -ipdb {refpdb} -O {ref}')

""
# load in reference molecule 
mol = Chem.MolFromMol2File(f'{ref}')
ref_feat = MACCSkeys.GenMACCSKeys(mol)

# load in ligand library
suppl = Chem.SDMolSupplier(f'{database}')

ms = [lig for lig in suppl if lig is not None]
lig_name, score = [], []

for lig in ms:
    Chem.AddHs(lig)
    name = lig.GetProp("_Name")
    print("Now processing ", name)
    lig_name.append(name)
    
    # generate pharmacophore features 
    feat = MACCSkeys.GenMACCSKeys(lig)
    
    # calculate TANIMOTO similarity
    s = DataStructs.FingerprintSimilarity(ref_feat,feat) # TANIMOTO
    score.append('{:.4f}'.format(s))

ref_name = [mol.GetProp("_Name").split('/')[-1]] * len(score)

rows = zip(ref_name,lig_name,score)
with open(f'{filename}.csv', "w") as f:
    writer = csv.writer(f, delimiter='\t')
    for row in rows:
        writer.writerow(row)
