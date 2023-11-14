import oddt
from oddt import toolkit
from oddt import docking
import glob
import os
import sys

lpath = sys.argv[1]

for sdf_path in glob.glob(f'{lpath}/*.sdf'):
   name = os.path.splitext(os.path.basename(sdf_path))[0]
  
   ligands = next(oddt.toolkit.readfile('sdf',sdf_path))
   oddt.docking.AutodockVina.write_vina_pdbqt(ligands,lpath,name_id=f'{str(name)}')
   os.remove(sdf_path)   