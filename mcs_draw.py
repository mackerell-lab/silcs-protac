from rdkit import Chem
from rdkit.Chem import Draw, rdFMCS, AllChem
import sys

# Load the molecules from input files
warhead = Chem.MolFromMolFile(sys.argv[1], sanitize=True, removeHs=True)         # removing Hs required for performance in certain cases
protac = Chem.MolFromMolFile(sys.argv[2], sanitize=True, removeHs=False) # leaving Hs required to minimize potentially confusing error

def draw_highlighted_warhead_in_protac(protac, warhead):
  # Find MCS between PROTAC and warhead
  mcs = rdFMCS.FindMCS([protac, warhead])
  mcs_mol = Chem.MolFromSmarts(mcs.smartsString)

  # Get matching atom indices in both molecules
  protac_match = protac.GetSubstructMatch(mcs_mol)
  warhead_match = warhead.GetSubstructMatch(mcs_mol)
  # Prepare lists of matching bonds
  protac_bond_matches = [
      bond.GetIdx() for bond in protac.GetBonds()
      if bond.GetBeginAtomIdx() in protac_match and bond.GetEndAtomIdx() in protac_match
  ]

  # Prepare drawing options
  AllChem.Compute2DCoords(protac)
  AllChem.Compute2DCoords(warhead)
  drawer_options = Draw.rdMolDraw2D.MolDrawOptions()
  drawer_options.highlightColor = (1, 0, 0)  # Highlight in red

  # Generate images for PROTAC with highlighted atoms and standalone warhead
  img_protac = Draw.MolToImage(protac, highlightAtoms=protac_match, highlightBonds=protac_bond_matches, options=drawer_options, size=(900, 900))
  img_warhead = Draw.MolToImage(warhead, size=(900, 900))

  return img_protac, img_warhead

im,im2 = draw_highlighted_warhead_in_protac(protac, warhead)
im.save('highlighted_protac.png')
im2.save('warhead.png')
