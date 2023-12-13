from rdkit import Chem
from IPython.display import SVG
from rdkit.Chem import AllChem, Draw, Descriptors, PandasTools
import os
import sys

s_name = sys.argv[1]
smiles_data = sys.argv[2]
lpath = sys.argv[3]

mh = Chem.AddHs(Chem.MolFromSmiles(str(smiles_data)))
AllChem.EmbedMolecule(mh, AllChem.ETKDGv2()) 
writer = Chem.SDWriter(f'{lpath}/{str(s_name)}.sdf')
writer.write(mh)

