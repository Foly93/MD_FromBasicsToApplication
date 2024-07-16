import sys
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Fragments
from rdkit.Chem.GraphDescriptors import BertzCT

in_file = sys.argv[1] # unfiltered.txt
out_file = sys.argv[2] # filtered.txt

# Read input smiles
smiles = pd.read_csv(in_file,header=None,
dtype=str,sep=' ').to_numpy(dtype=str)[:,0]
if 'smiles' in smiles:
    smiles = np.delete(smiles,0)

# Filtering Atoms that are not defined in later steps of the analysis
smi = [m for m in smiles if 'Si' not in m and '+' not in m and 'B' not in m and '-' not in m and 'I' not in m]

# Molecule should have one COO-, 4 or less Hbonddonors and acceptor, 1 imidazol or benzene group and BerzCT complexity < 500
ms = [Chem.MolFromSmiles(m) for m in smi]
ms = [x for x in ms if Fragments.fr_COO(x) == 1]
ms = [x for x in ms if Chem.rdMolDescriptors.CalcNumHBD(x) < 4]
ms = [x for x in ms if Chem.rdMolDescriptors.CalcNumHBA(x) < 4]
ms = [x for x in ms if Fragments.fr_imidazole(x) == 1 or Fragments.fr_benzene(x)==1]
ms = [x for x in ms if BertzCT(x) < 500 ]
sub_smiles = [ Chem.rdmolfiles.MolToSmiles(x) for x in ms ]

# Output
np.savetxt(out_file,sub_smiles,fmt='%s')
