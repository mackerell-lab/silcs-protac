#!/home/enord/local/miniconda3/envs/silcsprotac/bin/python

from rdkit import Chem
from rdkit.Chem import AllChem, SDWriter, rdFMCS
import argparse
import os

# Argument parser for input files
parser = argparse.ArgumentParser(prog='mcs_align.py', description='This script performs Maximum Common Substructure Alignment of warhead to ref_ligands, using the MCS to assign the coordinates of ref_ligands to warhead, and falling back on RMSD rigid alignment if MCS Alignment fails.')
parser.add_argument('-w', '--warhead', required=True, help='Warhead or fragment common among ref_ligands, single ligand in sdf/mol format')
parser.add_argument('-r', '--ref_ligands',  required=True, help='Ref coordinates for alignment containing warhead, can be multiple ligands in sdf format')
parser.add_argument('-o', '--output', default='output.sdf', help='output filename; default = output.sdf')
args = parser.parse_args()

# Load the molecules from input files
warhead = Chem.MolFromMolFile(args.warhead, sanitize=True, removeHs=True)         # removing Hs required for performance in certain cases
ref_ligands = Chem.SDMolSupplier(args.ref_ligands, sanitize=True, removeHs=False) # leaving Hs required to minimize potentially confusing error

warhead_basename = os.path.basename(args.warhead[:-4])
reflig_basename  = os.path.basename(args.ref_ligands)
output_filename = args.output
writer = SDWriter(output_filename)

def mcs_align(mcsmol, reflig, refmatch):
  rwmol = Chem.RWMol(mcsmol) # create new molecule based on MCS
  rwconf = Chem.Conformer(rwmol.GetNumAtoms()) # initialize conformer
  matches = rwmol.GetSubstructMatch(mcsmol)    # get the substructure matches, most/all of them

  # now loop over the reference ligand atom positions rdkit style
  refconf = reflig.GetConformer() 
  for j, match in enumerate(matches):
    # in the initialized conformer, set the atom positions of index match to the reference coords of atom index refmatch[j]
    # two lists of matching atom indices; (new, warhead-like substructure) match and (reference) refmatch[j] are a pair
    rwconf.SetAtomPosition(match, refconf.GetAtomPosition(refmatch[j]))
  rwmol.AddConformer(rwconf)

  return rwmol
  
for i, reflig in enumerate(ref_ligands):
  if reflig is None: continue
  # Check for MCS matches and extract/write the fragments
  mcs = rdFMCS.FindMCS([warhead, reflig])
  mcs_mol = Chem.MolFromSmarts(mcs.smartsString)
  ref_match = reflig.GetSubstructMatch(mcs_mol)
  war_match = warhead.GetSubstructMatch(mcs_mol)
  mcs_frac = len(war_match)/warhead.GetNumHeavyAtoms()
    
  try:
    if ref_match and war_match:
      embed_ref = mcs_align(mcs_mol, reflig, ref_match) # has the warhead structure, but gets reference coords
      # these two lines below:
      mcs_embed_rms = AllChem.CalcRMS(embed_ref, reflig)
      rmsd_warhead = AllChem.CalcRMS(embed_ref, warhead)#, map=list(zip(ref_match, war_match)))

      warhead_save = embed_ref
      warhead_save.SetProp('MCS_RMSD_to_warhead', '%.2f' % rmsd_warhead)
      warhead_save.SetProp('MCS_EMBED_RMSD_expect_~0.0', '%.2f' % mcs_embed_rms)
      warhead_save.SetProp('MCS_FRAC_expect_~1.0', '%.2f' % mcs_frac)
      print('MCS_RMSD_to_Warhead = %.2f MCS_EMBED_RMSD = %.2f MCS_FRAC = %.2f'%(rmsd_warhead,mcs_embed_rms, mcs_frac))

    else:
      print('No valid MCS matches found in either molecule for entry %d.'%i)
  except ValueError:
    print('Something went wrong with constrained embedding, performing RMSD rigid alignment')
    rms = AllChem.AlignMol(warhead, reflig, atomMap=list(zip(war_match, ref_match)))
    warhead_save.SetProp('RIGID_RMSD_to_reference', '%.02f'%rms)

  warhead_save.SetProp('Entry', str(i))
  warhead_save.SetProp('_Name', reflig.GetProp('_Name'))
  warhead_save.SetProp('Warhead name', warhead.GetProp('_Name'))
  writer.write(warhead_save)
     
writer.close()
