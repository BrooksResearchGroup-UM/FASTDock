# pycharmm FASTDock test script 
import sys
import shutil
import os
import numpy as np

# find the current working directory to use below
cwd = os.getcwd()
sys.path.append(f'{cwd}/script')
from FASTDock import dock_probe

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

## Set FASTDock parameters
exec(open(f'{cwd}/script/parameters.py').read())

# Establish the directory for the test run
testdir = f'{cwd}/{testdir}'
if not os.path.exists(testdir): os.mkdir(testdir)
os.chdir(testdir)

""
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
# #                           ##
# # GENERATE GRID FOR DOCKING ##
# #                           ##
# ##############################

## Build protein
read.psf_card(protein_psf, append = True)
read.pdb(protein_pdb, resid = True)
coor.orient(by_rms=True,by_mass=False,by_noro=False)

write.coor_pdb('final_protein.pdb')

# set box size and center
stats = coor.stat()
xs = stats['xmax'] - stats['xmin']
ys = stats['ymax'] - stats['ymin']
zs = stats['zmax'] - stats['zmin']

smax = np.amax([xs,ys,zs]) + buff

# set grid parameters, nonbonded options and generate grid
nbond_opt = {'cutnb' : 14 , 'ctofnb' : 12 , 
                                        'ctonnb' : 10, 'eps' : 2,
                                        'rdie' : True, 'e14fac' : 1}
pycharmm.NonBondedScript(**nbond_opt).run()
genGrid = grid.CDOCKER()
settings = {'xCen' : stats['xave'], 'yCen' : stats['yave'], 'zCen' : stats['zave'],
            'xMax' : smax, 'yMax' : smax, 'zMax' : smax,
            'dGrid' : dgrid,
            'emax' : 2, 'maxe' : 40, 'mine' : -20, 'dielec' : 2, 'flag_gpu' : True,
            'rdie' : True,
            'nbond_opt' : nbond_opt,
            'gridFile' : 'grid/protein_grid.bin',
            'probeFile' : f'"{cwd}/data/toppar/fftdock_c36prot_cgenff_probes.txt"'}

genGrid.setVar(settings)

if os.path.exists('grid'): shutil.rmtree('grid')
os.mkdir('grid')

genGrid.generate()

psf.delete_atoms(pycharmm.SelectAtoms().all_atoms())

###############################
# #   DOCK PROBES W/ FFTDOCK  ##
# ##############################

if os.path.exists('scores'): shutil.rmtree("scores")
if os.path.exists('pdb_dock'): shutil.rmtree('pdb_dock')
os.mkdir('scores')
os.mkdir('pdb_dock')

dock_probe(cwd, addmin, randlig, fragnames, nsave, unit=30)
