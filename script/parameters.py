# Modifiable parameters for the FASTDock workflow 

# Define the directory in which test is to be run
testdir = 'test_jupyter'
# Specify the path to the pyCHARMM library (for setting CHARMM_LIB_DIR)
pyCHARMM_LIB = '/home/brookscl/charmm/c47-dev-release/install-testpycharmm/lib'
# Set FASTDock Parameters 
# Note we would generally save 2000 but for expedience in this example
# we use 100
nsave = 100 # number of saved probe poses
# The default for minsize would be 10, but for this test we use a smaller number
# since we are only saving 20 docked poses.
minsize = 10 # number of members in clusters to advance
# Generally, we would minimize in grid on cpu but for epedience in example
# we don't.
addmin = True
randlig = True # randomize quaternion positions
dgrid = 1 # grid spacing
buff = 10 # grid buffer around protein
segid = 'PRO0' # segment ID for protein

fragnames = [ 'etha', 'pro2', 'etoh', 'tboh', 'acn' , 'mam1', 'dmf' , 'dmee', 'bald', 'benz', 'chxe', 'phen', 'acem', 'aco' , 'aald', 'urea', 'acet', 'mamm']
nfrag = len(fragnames) # number of probes 

# Note, embedding file path string in double quotes allows CHARMM to correctly 
# translate mixed case file names.
protein_pdb = f'"{cwd}/data/pdb/protein.pdb"' # path to protein pdb file
protein_psf = f'"{cwd}/data/pdb/protein.psf"' # path to protein psf file
