#!/home/enord/local/miniconda3/envs/silcsprotac/bin/python

from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS
import argparse
import os

# Argument parser for input files
parser = argparse.ArgumentParser(prog='mcs_lgfe.py', description='This script obtains the Maximum Common Substructure a warhead within ref_ligands, then sums up and writes out the GFEs of each ref_ligand atom in the MCS, giving an effective warhead LGFE.')
parser.add_argument('-w', '--warhead', required=True, help='Warhead or fragment common among ref_ligands, single ligand in sdf/mol format')
parser.add_argument('-r', '--ref_ligands',  required=True, help='Ref coordinates for alignment containing warhead, can be multiple ligands in sdf format')
args = parser.parse_args()

# Load the molecules from input files
warhead = Chem.MolFromMolFile(args.warhead, removeHs=True, sanitize=True)
ref_ligands = Chem.SDMolSupplier(args.ref_ligands, removeHs=False, sanitize=True)

warhead_basename = os.path.basename(args.warhead[:-4])
reflig_basename  = os.path.basename(args.ref_ligands)

def mcs_lgfe(reflig, refmatches):
  gfe_string = reflig.GetProp('atom.dprop.GFE')
  gfe_array  = [float(gfe) for gfe in gfe_string.split()]
  lgfe=0
  for j, match in enumerate(refmatches):
    lgfe += gfe_array[match]
  return lgfe
  
for i, reflig in enumerate(ref_ligands):
  mcs       = rdFMCS.FindMCS([warhead, reflig])
  mcs_mol   = Chem.MolFromSmarts(mcs.smartsString)
  ref_match = reflig.GetSubstructMatch(mcs_mol)
  war_match = warhead.GetSubstructMatch(mcs_mol)
    
  lgfe = mcs_lgfe(reflig, ref_match) # has the warhead structure, but gets reference coords
  mcs_frac = len(war_match)/warhead.GetNumAtoms()
  print('LGFE = %.2f, MCS_FRAC = %.2f'%(lgfe, mcs_frac))
