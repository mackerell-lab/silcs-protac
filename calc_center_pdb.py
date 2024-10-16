#!/usr/bin/env python
import MDAnalysis as mda
import numpy as np
import sys

u = mda.Universe(sys.argv[1], quiet=True)
lig = u.select_atoms('not name H')

com = lig.center_of_mass()
print(f'{com[0]:.1f},{com[1]:.1f},{com[2]:.1f}')
