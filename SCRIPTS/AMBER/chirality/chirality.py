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

def chiral_candidates(atoms):
  """
  Returns those atoms which have 4 neighbours in the atom graph (and thus might be chiral).

  :param atoms: atom graph
  :return: list of read_amber.Atom objects which have 4 neighbours in the atoms graph
  """
  candidates = [atom for atom in atoms.nodes() if len(nx.edges(atoms, atom)) == 4]
  return candidates

def not_in_dict(atoms):
  """
  Writes an error message to `chirality_error` and stderr, if an atom in the topology file is not present in the
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
    with open("chirality_warning", "w") as error_file:
      error_file.write("The following atoms were not present in the defined_valences dictionary\n")
      error_file.write("of the chirality script, and given a valence of 0 (not bonded):\n\n")
      sys.stderr.write("The following atoms were not present in the defined_valences dictionary\n")
      sys.stderr.write("of the chirality script, and given a valence of 0 (not bonded):\n\n")
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

def chiral_order(atoms, chiral_atom, depth=6):
  # Create a list of ordered atoms to be passed back
  ordered = []
  # Do a quick check whether there are multiple hydrogens
  neighbors = atoms.neighbors(chiral_atom)
  hydrogens = [atom for atom in neighbors if atom.element == "H"]
  if len(hydrogens) < 2:
    tree = nx.bfs_tree(atoms, chiral_atom)
    # Generate the list of shortest paths in the molecule, neglecting the trivial path [chiral_atom]
    paths = sorted(nx.single_source_shortest_path(tree, chiral_atom, depth).values(), reverse = True)[:-1]
    while paths:
      # Pop the first element (highest priority path) from the list of paths and remove any duplicates.
      path = paths.pop(0)
      paths_no_dups = [unpruned for unpruned in paths if unpruned != path]
      # If there are any duplicates, the paths list will be smaller and we can't resolve a highest priority yet.
      if len(paths_no_dups) != len(paths):
        paths = paths_no_dups
      # Otherwise, the path is higher priority than all the other paths, so its second atom is the neighbour with
      # highest priority.
      else:
        ranked_atom = path[1]
        ordered.append(ranked_atom)
        # Drop all the paths containing our ranked atom.
        paths = [unpruned for unpruned in paths if unpruned[1] is not ranked_atom]
  return ordered

def write_chirality_file(input_filename, output_filename, human_readable_filename):
  molecule = ra.parse_topology_file(input_filename)
  atoms = molecule.atoms
  chiral_cands = chiral_candidates(atoms)
  multi_bonds(atoms)
  chiral_centres = {}
  for i, chiral_atom in enumerate(chiral_cands):
    ordered = chiral_order(atoms, chiral_atom)
    if len(ordered) == 4:
      chiral_centres[chiral_atom] = ordered
  with open(output_filename, "w") as output_file:
    for atom in sorted(chiral_centres.keys(), cmp=lambda x, y: cmp(x.index, y.index)):
      # Write out the list of chiral atoms and their CIP-ranked neighbours.
      output_string = "{0:>8d}{1:>8d}{2:>8d}{3:>8d}{4:>8d}\n".format(atom.index + 1, 
                                                                     *[other_atom.index + 1 for other_atom in chiral_centres[atom]])
      output_file.write(output_string)
  with open(human_readable_filename, "w") as human_file:
    human_file.write("Atom names given below are in the format: <index (from 0)> <element> <atom name>\n")
    human_file.write("e.g.: 12 C CA is the 13th atom, is a carbon and is named \"CA\"\n\n")
    output_string = "{0:^16s}{1:^16s}{2:^16s}{3:^16s}{4:^16s}\n".format("central atom", "top ranked", "2nd ranked", "3rd ranked", "lowest ranked")
    human_file.write(output_string)
    output_string = "{0:^16s}{1:^16s}{2:^16s}{3:^16s}{4:^16s}\n".format("============", "==========", "==========", "==========", "=============")
    human_file.write(output_string)
    for atom in sorted(chiral_centres.keys(), cmp=lambda x, y: cmp(x.index, y.index)):
      # Write out the list of chiral atoms and their CIP-ranked neighbours - but this time readable for humans.
      output_string = "{0:^16s}{1:^16s}{2:^16s}{3:^16s}{4:^16s}\n".format(str(atom), *[str(other_atom) for other_atom in chiral_centres[atom]])
      human_file.write(output_string)

if __name__ == "__main__":
  write_chirality_file(sys.argv[1], ".chirality_list", "chirality_readable")
