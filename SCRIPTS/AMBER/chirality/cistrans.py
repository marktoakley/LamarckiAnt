#!/usr/bin/python

import read_amber as ra
import networkx as nx
import sys
from collections import defaultdict

class GhostAtom(ra.Atom):
  pass

# Set the default valence to 0, if the element isn't present in the defined_valences dictionary.
valences = defaultdict(int)
defined_valences = { "H": 1,
             "C": 4,
             "N": 3,
             "O": 2,
             "S": 2,
             "Se": 2,
             "P": 5 }

for k, v in defined_valences.items():
    valences[k] = v

def peptide_bonds(atoms):
  peptide_bond_list = []
  amide_carbons = []
# Identify amide carbons.
  for carbon in [atom for atom in atoms if atom.element is "C"]:
    if isinstance(carbon, GhostAtom):
      continue
    neighbors = atoms.neighbors(carbon)
    nitrogens = [neighbor for neighbor in neighbors if neighbor.element is "N"]
    oxygens = [neighbor for neighbor in neighbors if neighbor.element is "O" and 
               not isinstance(neighbor, GhostAtom)]
    ghost_oxygens = [neighbor for neighbor in neighbors if neighbor.element is "O" and 
                     isinstance(neighbor, GhostAtom)]
    if len(nitrogens) == 1 and len(oxygens) == 1 and len(ghost_oxygens) == 1:
      amide_carbons.append(carbon)
# Identify the O - C - N - H atoms.
  for carbon in amide_carbons:
    neighbors = atoms.neighbors(carbon)
    oxygen = [neighbor for neighbor in neighbors if neighbor.element is "O"][0]
    nitrogen = [neighbor for neighbor in neighbors if neighbor.element is "N"][0]
    # Proline doesn't have N-H in its peptide bond
    if nitrogen.residue.name == 'PRO':
        continue
    n_neighbors = atoms.neighbors(nitrogen)
    hydrogens = [neighbor for neighbor in n_neighbors if neighbor.element is "H"]
    # Skip amide groups (where the N has two H neighbours)
    if len(hydrogens) > 1:
        continue
    hydrogen = hydrogens[0]
    peptide_bond_list.append((oxygen, carbon, nitrogen, hydrogen))
  return peptide_bond_list

def not_in_dict(atoms):
  """
  Writes an error message to `cistrans_warning` and stderr, if an atom in the topology file is not present in the
  `defined_valences` dictionary.

  :param atoms: atom graph
  """
  missing = []
  for atom in atoms.nodes():
    try:
      _ = defined_valences[atom.element]
    except KeyError:
      missing.append(str(atom))
  if missing:
    with open("cistrans_warning", "w") as error_file:
      error_file.write("The following atoms were not present in the defined_valences dictionary\n")
      error_file.write("of the cistrans script, and given a valence of 0 (not bonded):\n\n")
      sys.stderr.write("The following atoms were not present in the defined_valences dictionary\n")
      sys.stderr.write("of the cistrans script, and given a valence of 0 (not bonded):\n\n")
      for atom_str in missing:
        error_file.write(atom_str + "\n")
      for atom_str in missing:
        sys.stderr.write(atom_str + "\n")

def multi_bonds(atoms):
  not_in_dict(atoms)
  multibonded = [atom for atom in atoms.nodes() 
                 if len(nx.edges(atoms, atom)) < valences[atom.element]]
  for i, atom in enumerate(multibonded):
    paired = False
    for other in atoms.neighbors(atom):
      if isinstance(other, GhostAtom):
        paired = True
        continue
      if len(nx.edges(atoms, other)) < valences[other.element]:
        ghost_atom = GhostAtom(**(atom.__dict__))
        ghost_atom.name = atom.name + "*"
        ghost_other = GhostAtom(**(other.__dict__))
        ghost_other.name = other.name + "*"
        atoms.add_edge(other, ghost_atom)
        atoms.add_edge(atom, ghost_other)
        paired = True
#    if not paired:
#      print "Could not find pair for atom", atom, atom.residue, atom.charge, nx.edges(atoms, atom)
  
def write_cis_trans_file(input_filename, output_filename):
  molecule = ra.parse_topology_file(input_filename)
  atoms = molecule.atoms
#  for atom in atoms:
#    print atom.residue
  multi_bonds(atoms)
  with open(output_filename, "w") as output_file:
    for bond in sorted(peptide_bonds(atoms), cmp=lambda x, y: cmp(x[0].index, y[0].index)):
    # Write out the list of atoms in peptide bonds (O - C - N - H).
      output_string = "{0:>8d}{1:>8d}{2:>8d}{3:>8d}\n".format(*map(lambda x: x.index + 1, bond))
      output_file.write(output_string)

if __name__ == "__main__":
  write_cis_trans_file(sys.argv[1], ".cis_trans_list")
