#!/home/enord/local/miniconda3/envs/silcsprotac/bin/python

from rdkit import Chem
from rdkit.Chem import AllChem, SDWriter, rdFMCS
import argparse
import os

# Argument parser for input files
parser = argparse.ArgumentParser(prog='mcs_align.py', description='This script performs Maximum Common Substructure Alignment of warhead to ref_ligands, using the MCS to assign the coordinates of ref_ligands to warhead, and falling back on RMSD rigid alignment if MCS Alignment fails.')
parser.add_argument('-w', '--warhead', required=True, help='Warhead or fragment common among ref_ligands, single ligand in sdf/mol format')
parser.add_argument('-r', '--ref_ligands',  required=True, help='Ref coordinates for alignment containing warhead, can be multiple ligands in sdf format')
parser.add_argument('-p', '--protonate', action='store_true', help='Protonate output coordinates (necessary for computing LGFE; note this will neutralize molecule)')
parser.add_argument('-o', '--output', default='output.sdf', help='output filename; default = output.sdf')
args = parser.parse_args()

# Load the molecules from input files
warhead = Chem.MolFromMolFile(args.warhead, removeHs=True, sanitize=True)
ref_ligands = Chem.SDMolSupplier(args.ref_ligands, removeHs=True, sanitize=True)
Chem.SanitizeMol(warhead, sanitizeOps=Chem.SanitizeFlags.SANITIZE_SETAROMATICITY)


warhead_basename = os.path.basename(args.warhead[:-4])
reflig_basename  = os.path.basename(args.ref_ligands)
output_filename = args.output
writer = SDWriter(output_filename)

def protonate_constrained_embed(mol, ref):
  molH = Chem.AddHs(mol)
  AllChem.MMFFSanitizeMolecule(molH)#, mmffVariant='MMFF94')
  AllChem.MMFFOptimizeMolecule(molH, mmffVariant='MMFF94')#, maxIters=1000, nonBondedThresh=50.0)
  return AllChem.ConstrainedEmbed(molH, ref)

def mcs_align(mcsmol, reflig, refmatch):
  rwmol = Chem.RWMol(mcsmol)
  rwconf = Chem.Conformer(rwmol.GetNumAtoms())
  matches = rwmol.GetSubstructMatch(mcsmol)

  refconf = reflig.GetConformer()
  for j, match in enumerate(matches):
    # Added atom position information from reference molecule
    rwconf.SetAtomPosition(match, refconf.GetAtomPosition(refmatch[j]))
  rwmol.AddConformer(rwconf)

  return rwmol
  
for i, reflig in enumerate(ref_ligands):
  if reflig is None: continue
  Chem.SanitizeMol(reflig, sanitizeOps=Chem.SanitizeFlags.SANITIZE_SETAROMATICITY)
  # Check for MCS matches and extract/write the fragments
  mcs = rdFMCS.FindMCS([warhead, reflig])
  mcs_mol = Chem.MolFromSmarts(mcs.smartsString)
  ref_match = reflig.GetSubstructMatch(mcs_mol)
  war_match = warhead.GetSubstructMatch(mcs_mol)
    
  try:
    if ref_match and war_match:
      embed_ref = mcs_align(mcs_mol, reflig, ref_match) # has the warhead structure, but gets reference coords
      # these two lines below:
      mcs_rms = AllChem.CalcRMS(embed_ref, reflig)#, map=list(zip(war_match, ref_match)))

      if args.protonate: warhead = protonate_constrained_embed(warhead, embed_ref)
      else: warhead = embed_ref
      mcs_frac = len(war_match)/warhead.GetNumHeavyAtoms()
      warhead.SetProp('MCS_RMSD_to_reference', '%.2f' % mcs_rms)
      warhead.SetProp('MCS_FRAC_to_reference', '%.2f' % mcs_frac)
      print('RMSD = %.2f, frac. atoms in MCS = %.2f'%(mcs_rms, mcs_frac))

    else:
      print('No valid MCS matches found in either molecule for entry %d.'%i)
  except ValueError:
    print('Something went wrong with constrained embedding, performing RMSD rigid alignment')
    rms = AllChem.AlignMol(warhead, reflig, atomMap=list(zip(war_match, ref_match)))
    if args.protonate: warhead = protonate_constrained_embed(warhead, warhead)
    warhead.SetProp('RIGID_RMSD_to_reference', '%.02f'%rms)

  warhead.SetProp('Entry', str(i))
  warhead.SetProp('_Name', reflig.GetProp('_Name'))
  warhead.SetProp('Warhead name', warhead.GetProp('_Name'))
  writer.write(warhead)
  print('Warhead ligand aligned to protac entry %d and written to %s'%(i,output_filename))
     
writer.close()
if args.protonate: print('WARNING: protonate mode results in worse embeddings and unstable results, proceed with extra caution!')
