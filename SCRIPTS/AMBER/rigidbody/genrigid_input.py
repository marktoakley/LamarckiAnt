# SCRIPT TO DEFINE THE RIGIDIFICATION FOR PROTEINS
#by Konstantin Roeder
#creates input files rbodyconfig and coordsinirigid and additionally a log file rigid.log

#! /usr/bin/env python

import sys
import time
import math

if len(sys.argv) > 1:
    print 'Script creates input files for RIGIDINIT keyword, needs inpcrd and pdb'
    print
    pdb_inp = sys.argv[1]
    coords_inp = sys.argv[2]
    #defines the tolerance for the checks whether rigid bodies are linear
    try:
        tolerance = float(sys.argv[3])
    except IndexError:
        tolerance = 0.01
    #Pymol use; change to False if pymol is not installed/used
    try:
        if sys.argv[4] == 'pymol':
            pymol_check = True
    except IndexError:
        pymol_check = False

else:
    tolerance = 0.01
    pymol_check = False
    pdb_inp = raw_input('PDB file: ')
    coords_inp = raw_input('Coords input: ')

if pymol_check:
    import __main__

    __main__.pymol_argv = ['pymol', '-qix']  # Pymol: quiet and no GUI(internal and external)
    import pymol

#class containing the main methods for the script, new functionality should be embedded here
class protein():
    def __init__(self, atom_dic, res_dic):
        self.atom_dic = atom_dic  #dictionary of all atoms with x,y,z coordinates and res id and name and atom name
        self.res_dic = res_dic    #dictionary containing a list of atom ids for all residues

    def num_atoms(self): #number of atoms
        return len(self.atom_dic)

    def num_res(self): #number of residues
        return len(self.res_dic)

    #returns the name of a given residue
    def get_residue_name(self, residue):
        return atom_dic[res_dic[residue][0]][1]

    #prints the atoms in a specified residue(need number of residue)
    #includes the atom number, atom name and coordinates
    def get_atoms_in_res(self, residue):
        all_atoms = self.atom_dic.keys()
        atoms_in_res = []
        for i in range(1, len(self.atom_dic) + 1):
            atom = self.atom_dic[i]
            if residue == atom[2]:
                atoms_in_res.append(all_atoms[i - 1])
                residuename = atom[1]
        if atoms_in_res == []:
            return 'Residue not in molecule'
        print 'Residue:', residue, residuename
        print
        print 'Atoms:'
        for i in range(0, len(atoms_in_res)):
            j = atoms_in_res[i]
            atom = self.atom_dic[j]
            print j, atom[0], atom[3], atom[4], atom[5]
        return

    #returns a list of atoms in a residue
    def get_atom_list(self, residue):
        all_atoms = self.atom_dic.keys()
        atomlist = []
        for i in range(1, len(self.atom_dic) + 1):
            atom = self.atom_dic[i]
            if residue == atom[2]:
                atomlist.append(all_atoms[i - 1])
        return atomlist

    #allows the input of a centre from the user, is complemented by the selection of the CoM as centre
    def choose_centre(self):
        check = raw_input('Display atoms in residue (y/n)?  ')
        if check == 'y':
            while 1:  #acts as input mask --> once entered any number of residues can be viewed at
                residue_number = raw_input('Residue number: (enter x to leave) ')
                if residue_number == 'x':
                    print
                    break
                else:
                    try:
                        self.get_atoms_in_res(int(residue_number))
                        print
                    except ValueError:
                        pass
        while 1:  # loop for choice of molecule, only accepts atoms in the molecule
            centre_input = int(raw_input('Choose atom:  '))
            if centre_input > self.num_atoms() or centre_input < 1:
                print 'Not in molecule.'
            else:
                break
        return centre_input

    #transforms a given list of atoms into a list of residues
    #can be used to find the residues in a sphere from the output of the atoms in sphere method
    def transform_atomlist_to_reslist(self, atomlist):
        reslist = []
        for i in atomlist:
            atom_info = self.atom_dic[i]
            k = 0
            if reslist != []:
                for j in range(0, len(reslist)):
                    if reslist[j] == atom_info[2]:
                        k = 1
            if k == 0:
                reslist.append(atom_info[2])
            else:
                pass
        return reslist

    def transform_reslist_to_atomlist(self, reslist):
        atomlist = []
        for i in reslist:
            atomlist += self.res_dic[i]
        return atomlist

    #get residue name from id
    def res_name_from_id(self, residue_id):
        atom_id = res_dic[residue_id][0]
        return atom_dic[atom_id][1]

    #checks for unknown residues/ligands, add standard residues here if needed
    def get_ligands_unknown_res(self):
        list_res_names = ['ALA', 'ARG', 'ASN', 'ASP', 'ASX', 'CYS', 'CYX', 'GLN', 'GLU', 'GLY', 'GLX', 'HIS', 'HIE',
                          'HID', 'HIP', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
                          'NALA', 'NARG', 'NASN', 'NASP', 'NASX', 'NCYS', 'NCYX', 'NGLN', 'NGLU', 'NGLY', 'NGLX',
                          'NHIS', 'NHIE', 'NHID', 'NHIP', 'NILE', 'NLEU', 'NLYS', 'NMET', 'NPHE', 'NPRO', 'NSER',
                          'NTHR', 'NTRP', 'NTYR', 'NVAL', 'CALA', 'CARG', 'CASN', 'CASP', 'CASX', 'CCYS', 'CCYX',
                          'CGLN', 'CGLU', 'CGLY', 'CGLX', 'CHIS', 'CHIE', 'CHID', 'CHIP', 'CILE', 'CLEU', 'CLYS',
                          'CMET', 'CPHE', 'CPRO', 'CSER', 'CTHR', 'CTRP', 'CTYR', 'CVAL']
        ligand_dic = {}
        for i in self.res_dic.keys():
            if self.res_name_from_id(i) not in list_res_names:
                ligand_dic[i] = self.res_dic[i]
        return ligand_dic

    #returns a list of atoms within a sphere of entered radius
    #needs an atom number as centre
    def get_atoms_in_sphere(self, centre_atom, radius):
        all_atoms = self.atom_dic.keys()
        atomlist = []
        for i in range(1, len(self.atom_dic) + 1):
            atom = self.atom_dic[i]
            distance = (atom[3] - centre_atom[0]) ** 2 + (atom[4] - centre_atom[1]) ** 2 + (atom[5] - centre_atom[
                2]) ** 2
            if distance <= radius:
                atomlist.append(all_atoms[i - 1])
            else:
                pass
        return atomlist

    #returns the atoms within an ellipsoidal volume
    def get_atoms_in_ellipsoid(self, centre_atom, a, b, c):
        all_atoms = self.atom_dic.keys()
        atomlist = []
        for i in range(1, len(self.atom_dic) + 1):
            atom = self.atom_dic[i]
            x1 = ((atom[3] - centre_atom[0]) / a) ** 2
            x2 = ((atom[4] - centre_atom[1]) / b) ** 2
            x3 = ((atom[5] - centre_atom[2]) / c) ** 2
            distance = x1 + x2 + x3
            if distance <= 1:
                atomlist.append(all_atoms[i - 1])
            else:
                pass
        return atomlist

    #returns a list of atoms in a residue
    def get_atom_list(self, residue):
        all_atoms = self.atom_dic.keys()
        atomlist = []
        for i in range(1, len(self.atom_dic) + 1):
            atom = self.atom_dic[i]
            if residue == atom[2]:
                atomlist.append(all_atoms[i - 1])
        return atomlist

    #gives the mass of any residue containing only C, N, O, H, D and S, other masses can only be used if they are added to the dictionary mass_dic below
    def res_mass(self, residue):
        atom_list = self.get_atom_list(residue)
        mass_dic = {'H': 1.007825, 'D': 2.014102, 'C': 12.0116, 'N': 14.00728, 'O': 15.99977, 'S': 32.076}
        mass = 0
        for i in atom_list:
            for j in mass_dic.keys():
                if j == self.atom_dic[i][0][0]:
                    mass += mass_dic[j]
        return mass

    #gives the coordinates for the mass_weighted centre of a residue
    def mass_weighted_centre(self, residue):
        atom_list = self.get_atom_list(residue)
        mass_dic = {'H': 1.007825, 'D': 2.014102, 'C': 12.0116, 'N': 14.00728, 'O': 15.99977, 'S': 32.076}
        X = 0
        Y = 0
        Z = 0
        mass_res = float(self.res_mass(residue))
        for i in atom_list:
            for j in mass_dic.keys():
                if j == self.atom_dic[i][0][0]:
                    X += mass_dic[j] * self.atom_dic[i][3]
        for i in atom_list:
            for j in mass_dic.keys():
                if j == self.atom_dic[i][0][0]:
                    Y += mass_dic[j] * self.atom_dic[i][4]
        for i in atom_list:
            for j in mass_dic.keys():
                if j == self.atom_dic[i][0][0]:
                    Z += mass_dic[j] * self.atom_dic[i][5]
        X = X / mass_res
        Y = Y / mass_res
        Z = Z / mass_res
        print 'Mass weighted centre: ', (X, Y, Z)
        return (X, Y, Z)


    #returns all residues of a list of type 'name'
    def get_all_res_name(self, res_input, name):
        res_list = []
        for i in res_input:
            if name == self.get_residue_name(i):
                res_list.append(i)
        return res_list

    #atom id from atom name given a residue id
    def get_atom_id_from_name(self, atom_name, res_num):
        for i in self.res_dic[res_num]:
            if self.atom_dic[i][0] == atom_name:
                return i

    #dot product of two vectors
    def dot_product(self, vector1, vector2):
        return (vector1[0] * vector2[0] + vector1[1] * vector2[1] + vector1[2] * vector2[2])

    #magnitude of a vector
    def magnitude(self, vector):
        return math.sqrt(vector[0] ** 2 + vector[1] ** 2 + vector[2] ** 2)

    #angle between three atoms
    def get_angle(self, (atom1, atom2, atom3)):
        vector1 = (self.atom_dic[atom2][3] - self.atom_dic[atom1][3], self.atom_dic[atom2][4] - self.atom_dic[atom1][4]
                   , self.atom_dic[atom2][5] - self.atom_dic[atom1][5])
        vector2 = (self.atom_dic[atom3][3] - self.atom_dic[atom1][3], self.atom_dic[atom3][4] - self.atom_dic[atom1][4]
                   , self.atom_dic[atom3][5] - self.atom_dic[atom1][5])
        dot_product = self.dot_product(vector1, vector2)
        magnitudes = self.magnitude(vector1) * self.magnitude(vector2)
        angle = math.acos(dot_product / magnitudes)
        angle = math.degrees(angle)
        if (angle > 0.00 and angle <= 180.00):
            return angle
        else:
            return angle - 180.00

    #checks a list of atoms whether it is linear, creates a triple list of all possible combinations and then computes angles
    def check_linear(self, atom_list, tolerance):
        triple_list = []  # list of atom triples to check for linearity
        for atom1 in atom_list:
            for atom2 in atom_list:
                for atom3 in atom_list:
                    if (atom1 != atom2) and (atom1 != atom3) and (atom2 != atom3):
                        triple_list.append((atom1, atom2, atom3))
                    else:
                        pass
        for triple in triple_list:
            if (self.get_angle(triple) > (0.00 + tolerance)) and (self.get_angle(triple) < (180.00 - tolerance)):
                return False
        return True

    #returns a dictionary including all the information needed for the localised rigidification scheme,
    #if an extension is made to the default options for rigidification, this must change the local_scheme list!
    #the input scheme must be a list of 0,1,... for the possible options within the rigidification for a given residue
    #any additional scheme can be added below that will be available on default if chosen
    def get_rigid_for_local(self, res_list, local_scheme, atoms_used):
        #local_scheme=['PRO','ARG',('HIS','HIE','HID','HIP'),'LYS','ASP','ASN','GLU','GLN','PHE','TYR','TRP']
        groups_res = {}
        i = 1
        #all following cases are for a single residue name with given patterns
        #when adding additional ones check, that atoms can only be chosen once as no atom can be part of more than one rigid body
        if local_scheme[0] == 1:
            l_res = self.get_all_res_name(res_list, 'PRO')
            for j in l_res:
                atom_list = []
                for k in ['C', 'O', 'N', 'CD', 'CG', 'CB', 'CA']:
                    atom_list.append(self.get_atom_id_from_name(k, int(j)))
                groups_res[i] = [j, 'PRO'] + atom_list
                atoms_used += atom_list
                i += 1
            l_res = []
        if local_scheme[0] == 2:
            l_res = self.get_all_res_name(res_list, 'PRO')
            for j in l_res:
                atom_list = []
                for k in ['N', 'CD', 'CG', 'CB', 'CA']:
                    atom_list.append(self.get_atom_id_from_name(k, int(j)))
                groups_res[i] = [j, 'PRO'] + atom_list
                atoms_used += atom_list
                i += 1
            l_res = []
        if local_scheme[1] == 1:
            l_res = self.get_all_res_name(res_list, 'ARG')
            for j in l_res:
                atom_list = []
                for k in ['NE', 'HE', 'CZ', 'NH1', 'HH11', 'HH12', 'NH2', 'HH21', 'HH22']:
                    atom_list.append(self.get_atom_id_from_name(k, int(j)))
                groups_res[i] = [j, 'ARG'] + atom_list
                atoms_used += atom_list
                i += 1
            l_res = []
        if local_scheme[2] == 1:
            l_res = self.get_all_res_name(res_list, 'HIS')
            for j in l_res:
                atom_list = []
                for k in ['CG', 'ND1', 'CE1', 'NE2', 'CD2']:
                    atom_list.append(self.get_atom_id_from_name(k, int(j)))
                groups_res[i] = [j, 'HIS'] + atom_list
                atoms_used += atom_list
                i += 1
            l_res = []
            l_res = self.get_all_res_name(res_list, 'HIE')
            for j in l_res:
                atom_list = []
                for k in ['CG', 'ND1', 'CE1', 'NE2', 'CD2']:
                    atom_list.append(self.get_atom_id_from_name(k, int(j)))
                groups_res[i] = [j, 'HIE'] + atom_list
                atoms_used += atom_list
                i += 1
            l_res = []
            l_res = self.get_all_res_name(res_list, 'HID')
            for j in l_res:
                atom_list = []
                for k in ['CG', 'ND1', 'CE1', 'NE2', 'CD2']:
                    atom_list.append(self.get_atom_id_from_name(k, int(j)))
                groups_res[i] = [j, 'HID'] + atom_list
                atoms_used += atom_list
                i += 1
            l_res = []
            l_res = self.get_all_res_name(res_list, 'HIP')
            for j in l_res:
                atom_list = []
                for k in ['CG', 'ND1', 'CE1', 'NE2', 'CD2']:
                    atom_list.append(self.get_atom_id_from_name(k, int(j)))
                groups_res[i] = [j, 'HIP'] + atom_list
                atoms_used += atom_list
                i += 1
            l_res = []
        if local_scheme[3] == 1:
            l_res = self.get_all_res_name(res_list, 'LYS')
            for j in l_res:
                atom_list = []
                for k in ['NZ', 'HZ1', 'HZ2', 'HZ3']:
                    atom_list.append(self.get_atom_id_from_name(k, int(j)))
                groups_res[i] = [j, 'LYS'] + atom_list
                atoms_used += atom_list
                i += 1
            l_res = []
        if local_scheme[4] == 1:
            l_res = self.get_all_res_name(res_list, 'ASP')
            for j in l_res:
                atom_list = []
                for k in ['CB', 'CG', 'OD1', 'OD2']:
                    atom_list.append(self.get_atom_id_from_name(k, int(j)))
                groups_res[i] = [j, 'ASP'] + atom_list
                atoms_used += atom_list
                i += 1
            l_res = []
        if local_scheme[5] == 1:
            l_res = self.get_all_res_name(res_list, 'ASN')
            for j in l_res:
                atom_list = []
                for k in ['CB', 'CG', 'OD1', 'ND2', 'HD21', 'HD22']:
                    atom_list.append(self.get_atom_id_from_name(k, int(j)))
                groups_res[i] = [j, 'ASN'] + atom_list
                atoms_used += atom_list
                i += 1
            l_res = []
        if local_scheme[6] == 1:
            l_res = self.get_all_res_name(res_list, 'GLU')
            for j in l_res:
                atom_list = []
                for k in ['CG', 'CD', 'OE1', 'OE2']:
                    atom_list.append(self.get_atom_id_from_name(k, int(j)))
                groups_res[i] = [j, 'GLU'] + atom_list
                atoms_used += atom_list
                i += 1
            l_res = []
        if local_scheme[7] == 1:
            l_res = self.get_all_res_name(res_list, 'GLN')
            for j in l_res:
                atom_list = []
                for k in ['CG', 'CD', 'OE1', 'NE2', 'HE21', 'HE22']:
                    atom_list.append(self.get_atom_id_from_name(k, int(j)))
                groups_res[i] = [j, 'GLN'] + atom_list
                atoms_used += atom_list
                i += 1
            l_res = []
        if local_scheme[8] == 1:
            l_res = self.get_all_res_name(res_list, 'PHE')
            for j in l_res:
                atom_list = []
                for k in ['CG', 'CD1', 'HD1', 'CE1', 'HE1', 'CZ', 'HZ', 'CE2', 'HE2', 'CD2', 'HD2']:
                    atom_list.append(self.get_atom_id_from_name(k, int(j)))
                groups_res[i] = [j, 'PHE'] + atom_list
                atoms_used += atom_list
                i += 1
            l_res = []
        if local_scheme[9] == 1:
            l_res = self.get_all_res_name(res_list, 'TYR')
            for j in l_res:
                atom_list = []
                for k in ['CG', 'CD1', 'HD1', 'CE1', 'HE1', 'CZ', 'OH', 'CE2', 'HE2', 'CD2', 'HD2']:
                    atom_list.append(self.get_atom_id_from_name(k, int(j)))
                groups_res[i] = [j, 'TYR'] + atom_list
                atoms_used += atom_list
                i += 1
            l_res = []
        if local_scheme[10] == 1:
            l_res = self.get_all_res_name(res_list, 'TRP')
            for j in l_res:
                atom_list = []
                for k in ['CG', 'CD1', 'HD1', 'NE1', 'HE1', 'CE2', 'CZ2', 'HZ2', 'CH2', 'HH2', 'CZ3', 'CE3', 'HE3',
                          'CD2']:
                    atom_list.append(self.get_atom_id_from_name(k, int(j)))
                groups_res[i] = [j, 'TRP'] + atom_list
                atoms_used += atom_list
                i += 1
            l_res = []
        if local_scheme[10] == 2:
            l_res = self.get_all_res_name(res_list, 'TRP')
            for j in l_res:
                atom_list1 = []
                atom_list2 = []
                for k in ['CG', 'CD1', 'HD1', 'NE1', 'HE1']:
                    atom_list1.append(self.get_atom_id_from_name(k, int(j)))
                for k in ['CE2', 'CZ2', 'HZ2', 'CH2', 'HH2', 'CZ3', 'CE3', 'HE3', 'CD2']:
                    atom_list2.append(self.get_atom_id_from_name(k, int(j)))
                groups_res[i] = [j, 'TRP'] + atom_list1
                groups_res[i + 1] = [j, 'TRP'] + atom_list2
                atoms_used += atom_list1
                atoms_used += atom_list2
                i += 2
            l_res = []
        return (groups_res, atoms_used)


#class to implement a live feed in pymol - not yet tested
class PymolView():
    #initialises pymol and starts it without GUI, loads a given file and sets it to a white cartoon representation
    def __init__(self, file_name):
        print 'Launching Pymol'
        pymol.finish_launching()
        pymol.cmd.load(file_name)
        print '%s loaded in Pymol' % file_name

        #show molecule as cartoon and hide lines
        pymol.cmd.show('cartoon')
        pymol.cmd.hide('lines')
        pymol.cmd.set('cartoon_color', 0)

    def set_colour_range(self, start, finish, colour):
        #setting residue range to any colour
        pymol.cmd.set('cartoon_color', colour, 'resi %s-%s' % (str(start), str(finish)))

    #reset all of the molecule to a white colour
    #def reset_colour(self,length):
    #    pymol.cmd.set('cartoon_color',0,'resi 1-%s'%str(length))

    def set_colour_res(self, res, colour):
        pymol.cmd.set('cartoon_color', colour, 'resi %s' % str(res))

    def __del__(self):
        print 'Quit pymol'
        pymol.cmd.quit()

    def apply_colour_scheme(self, atomistic_list, local_list, colour_atomistic, colour_local):
        #self.reset_colour()
        for i in atomistic_list:
            self.set_colour_res(i, colour_atomistic)
        for j in local_list:
            self.set_colour_res(j, colour_local)


###FUNCTION TO PARSE THE PDB AND COORDS.INPCRD FILE ###################################################################
def reading_pdb_inpcrd(pdb, inpcrd):
    dic = {}
    try:
        atom_list = []
        with open(pdb, 'r') as f_pdb:
            for line in f_pdb:
                if line.split()[0] == 'ATOM':
                    atom_list.append((line[12:16].strip(), line[17:20].strip(), int(line[22:26])))
        with open(inpcrd, 'r') as f_crd:
            for _ in xrange(2):
                next(f_crd)
            i = 1
            for line in f_crd:
                atom_line = line.split()
                try:
                    atom_list[2 * i - 2] = atom_list[2 * i - 2] + (
                        float(atom_line[0]), float(atom_line[1]), float(atom_line[2]))
                except IndexError:
                    pass
                try:
                    atom_list[2 * i - 1] = atom_list[2 * i - 1] + (
                        float(atom_line[3]), float(atom_line[4]), float(atom_line[5]))
                    i = i + 1
                except IndexError:
                    pass
        for i in range(1, len(atom_list) + 1):
            dic[i] = atom_list[i - 1]
        return dic
    except IOError:
        return dic


def create_res_dic(atom_dic):
    res_dic = {}
    res_list = []
    for i in atom_dic.keys():
        res_list.append(atom_dic[i][2])
    res_list = list(set(res_list))
    for j in res_list:
        res_dic[j] = []
    for k in atom_dic.keys():
        res_dic[atom_dic[k][2]].append(k)
    return res_dic


###FUNCTIONS HELPING TO MANAGE THE LISTS AND DICTIONARY ###############################################################
def merge_lists(list_a, list_b):
    return sorted(list(set(list_a + list_b)))


def remove_selection_from_list(list_res, selection):
    new_res_list = []
    for i in list_res:
        if i not in selection:
            new_res_list.append(i)
    return sorted(new_res_list)


def check_res_exists(res_list, selection):
    new_selection = []
    for i in selection:
        if i in res_list:
            new_selection.append(i)
    return new_selection


def merge_new_and_old_res_lists(frozen, local, atomistic, sel_local, sel_atomistic):
    new_sel_atomistic = check_res_exists(frozen, sel_atomistic)
    new_sel_atomistic_from_local = check_res_exists(local, sel_atomistic)
    frozen = remove_selection_from_list(frozen, new_sel_atomistic)
    atomistic = merge_lists(atomistic, new_sel_atomistic)
    atomistic = merge_lists(atomistic, new_sel_atomistic_from_local)
    local = remove_selection_from_list(local, new_sel_atomistic_from_local)
    new_sel_local = check_res_exists(frozen, sel_local)
    local = merge_lists(local, new_sel_local)
    frozen = remove_selection_from_list(frozen, new_sel_local)
    return (sorted(atomistic), sorted(local), sorted(frozen))


def removal_selection(frozen, local, atomistic, sel_local, sel_atomistic):
    frozen = merge_lists(frozen, sel_local)
    frozen = merge_lists(frozen, sel_atomistic)
    local = remove_selection_from_list(local, sel_local)
    atomistic = remove_selection_from_list(atomistic, sel_atomistic)
    return (sorted(atomistic), sorted(local), sorted(frozen))


###FUNCTION TO SELECT RESIDUES ########################################################################################
def choose_res(mol):
    print '1 - choose sphere/ellipsoid of unfrozen atoms and add localised rigidification'
    print '2 - choose a range of residues or single residues to unfreeze or locally rigidify'
    print
    log = []
    try:
        method = int(raw_input('Method:  '))
    except IOError:
        method = 0
    unfrozen_res = []
    local_res = []

    if method == 0:
        print 'Invalid input'

    # ellipsoidal or spherical volume
    elif method == 1:
        while 1:  #choose the geometry: currently ellipsoidal 'e' or spherical 's'
            log.append('Method: sphere / ellipsoid')
            while 1:
                centre_method = raw_input('Mass weighted centre of residue (1) or atom in a selected residue (2)?  ')
                try:
                    if int(centre_method) == 1:
                        residue = raw_input('Residue:  ')
                        centre_xyz = mol.mass_weighted_centre(int(residue))
                        log.append('Mass-weighted centre: Residue' + residue + ' ' + mol.get_residue_name(
                            int(residue)) + '  ' + str(centre_xyz[0]) + '  ' + str(centre_xyz[1]) + '  ' + str(
                            centre_xyz[2]))
                        break
                    elif int(centre_method) == 2:
                        centre = mol.choose_centre()
                        centre_xyz = (mol.atom_dic[centre][3], mol.atom_dic[centre][4], mol.atom_dic[centre][5])
                        log.append('Centre: Atom  ' + str(centre) + '  ' + str(centre_xyz[0]) + '  ' + str(
                            centre_xyz[1]) + '  ' + str(centre_xyz[2]))
                        break
                    else:
                        print 'Invalid choice'
                except ValueError:
                    print 'Invalid input'

            radius_input = raw_input(
                'Radius of atomistic region (enter three values for an ellipsoid and one for a sphere):  ')
            l_radii = radius_input.split()
            try:
                if len(l_radii) == 1:
                    log.append('Radius of atomistic region: ' + l_radii[0])
                    l_atom = mol.get_atoms_in_sphere(centre_xyz, float(l_radii[0]))
                    l_res = mol.transform_atomlist_to_reslist(l_atom)
                    unfrozen_res += l_res
                    break
                elif len(l_radii) == 3:
                    log.append('Radii for the atomistic region: ' + l_radii[0] + '  ' + l_radii[1] + '  ' + l_radii[2])
                    l_atom = mol.get_atoms_in_ellipsoid(centre_xyz, float(l_radii[0]), float(l_radii[1]),
                                                        float(l_radii[2]))
                    l_res = mol.transform_atomlist_to_reslist(l_atom)
                    unfrozen_res += l_res
                    break
                else:
                    print 'Wrong input'
            except ValueError:
                print 'Invalid input'
        print
        while 1:  #geometry is same as for atomistic region, if a change is wanted employ the above
            check = raw_input('Add localised rigidification (y/n)?  ')
            if check == 'yes' or check == 'y':
                log.append('Local rigidification added')
                radius_input = raw_input(
                    'Radius of atomistic region (enter three values for an ellipsoid and one for a sphere):  ')
                l_radii_outer = radius_input.split()
                try:
                    if len(l_radii_outer) == 1:
                        log.append('Outer radius: ' + l_radii_outer[0])
                        l_atom_out = mol.get_atoms_in_sphere(centre_xyz, float(l_radii_outer[0]))
                        l_res_out = mol.transform_atomlist_to_reslist(l_atom_out)
                        for i in l_res_out:
                            if i not in l_res:
                                local_res.append(i)
                        break
                    elif len(l_radii_outer) == 3:
                        log.append(
                            'Radii for locally rigid region: ' + l_radii_outer[0] + '  ' + l_radii_outer[1] + '  ' +
                            l_radii_outer[2])
                        l_atom_out = mol.get_atoms_in_ellipsoid(centre_xyz, float(l_radii_outer[0]),
                                                                float(l_radii_outer[1]),
                                                                float(l_radii_outer[2]))
                        l_res_out = mol.transform_atomlist_to_reslist(l_atom_out)
                        for i in l_res_out:
                            if i not in l_res:
                                local_res.append(i)
                        break
                    else:
                        print 'Wrong input'
                except ValueError:
                    print 'Invalid input'
            else:
                break

    #rigidification using a range of residues, sanity check in function that chooses
    elif method == 2:
        log.append('Method:  range selection')
        check = raw_input('Define atomistic residues (y/n)?  ')
        if check == 'yes' or check == 'y':
            res = raw_input(
                'Enter residues separated by spaces, use r followed by two residues to define a range: ')  #input allows to enter 'r num1 num2' being the range between num1 and num2
            res = res.split()
            i = 0
            log_res = ''
            while i < len(res):
                if res[i] == 'r':  #if a range is indicated
                    i += 1
                    start_range = int(res[i])  #finds starting residue
                    i += 1
                    end_range = int(res[i])  #finds final residue
                    for j in range(start_range, end_range + 1):
                        unfrozen_res.append(j)
                        log_res += str(j) + ' '
                    i += 1
                else:
                    try:
                        unfrozen_res.append(int(res[i]))
                        log_res += res[i] + ' '
                        i += 1
                    except ValueError:
                        i += 1
            log.append('Residues entered as atomistically: ' + log_res)
        check = raw_input('Locally rigidify (y/n)?  ')
        if check == 'yes' or check == 'y':
            res = raw_input(
                'Enter residues separated by spaces, use r followed by two residues to define a range: ')  #input allows to enter 'r num1 num2' being the range between num1 and num2
            res = res.split()
            i = 0
            log_res = ''
            while i < len(res):
                if res[i] == 'r':  #if a range is indicated
                    i += 1
                    start_range = int(res[i])  #finds starting residue
                    i += 1
                    end_range = int(res[i])  #finds final residue
                    i += 1
                    for j in range(start_range, end_range + 1):
                        if j <= mol.num_res():
                            local_res.append(j)
                            log_res += str(j) + ' '
                        else:
                            break
                else:
                    try:
                        if int(res[i]) <= mol.num_res():
                            local_res.append(int(res[i]))
                            log_res += res[i] + ' '
                            i += 1
                    except ValueError:
                        i += 1
            log.append('Residues entered for local rigidification: ' + log_res)

    else:
        return 'Wrong input'
    return (unfrozen_res, local_res, log)


def remove_res(mol):
    unfrozen_res = []
    local_res = []
    log = []
    log.append('Removal of selected residues')
    check = raw_input('Remove from the atomistic residues (y/n)?  ')
    if check == 'yes' or check == 'y':
        res = raw_input(
            'Enter residues separated by spaces, use r followed by two residues to define a range: ')  #input allows to enter 'r num1 num2' being the range between num1 and num2
        res = res.split()
        i = 0
        log_res = ''
        while i < len(res):
            if res[i] == 'r':  #if a range is indicated
                i += 1
                start_range = int(res[i])  #finds starting residue
                i += 1
                end_range = int(res[i])  #finds final residue
                for j in range(start_range, end_range + 1):
                    unfrozen_res.append(j)
                    log_res += str(j) + ' '
                i += 1
            else:
                try:
                    unfrozen_res.append(int(res[i]))
                    log_res += res[i] + ' '
                    i += 1
                except ValueError:
                    i += 1
        log.append('Residues entered to be removed: ' + log_res + "\n")
    check = raw_input('Remove from the local selection (y/n)?  ')
    if check == 'yes' or check == 'y':
        res = raw_input(
            'Enter residues separated by spaces, use r followed by two residues to define a range: ')  #input allows to enter 'r num1 num2' being the range between num1 and num2
        res = res.split()
        i = 0
        log_res = ''
        while i < len(res):
            if res[i] == 'r':  #if a range is indicated
                i += 1
                start_range = int(res[i])  #finds starting residue
                i += 1
                end_range = int(res[i])  #finds final residue
                i += 1
                for j in range(start_range, end_range + 1):
                    if j <= mol.num_res():
                        local_res.append(j)
                        log_res += str(j) + ' '
                    else:
                        break
            else:
                try:
                    if int(res[i]) <= mol.num_res():
                        local_res.append(int(res[i]))
                        log_res += res[i] + ' '
                        i += 1
                except ValueError:
                    i += 1
        log.append('Residues entered to be removed from local rigidification: ' + log_res + "\n")
    return (unfrozen_res, local_res, log)


###VARIABLE DEFINITION ################################################################################################
atom_dic = {}  #dictionary containing all information sorted by atom number:
               #the information is in a tuple: (atomname, residue name, residue number, x, y, z)
res_dic = {}  #dictionary containing all atoms per residue sorted by residue number
res_atomistic = []  #residue list that are treated fully atomistic
res_local = []  #residue list that are treated locally rigid
res_frozen = []  #residue list that are treated fully rigid
atoms_used = []  #atoms used in rigid bodies, prevents double usage!
groups = {}  #dictionary of rigid groups
groups_local = {}
num_rbody = 0  #number of rigid bodies defined

###PARSING OF INPUT FILES #############################################################################################
while 1:  # starts input selection
    atom_dic = reading_pdb_inpcrd(pdb_inp, coords_inp)  # opens files and loads them into a dictionary
    res_dic = create_res_dic(atom_dic)
    if atom_dic != {}:  # if the dictionary contains data the data is read into the molecule class
        loaded_mol = protein(atom_dic, res_dic)
        print
        print 'Molecule was loaded'
        print 'Number of residues: ', loaded_mol.num_res()
        print 'Number of atoms: ', loaded_mol.num_atoms()
        print
        break
    else:
        print 'Files could not be found or they are empty - please try again'
        print
        pdb_inp = raw_input('PDB file: ')
        coords_inp = raw_input('Coords input: ')

check_file = open('rigid.log', 'w')
check_file.write('Input file for structure: ' + pdb_inp + "\n")
check_file.write('Coordinates file used: ' + coords_inp + "\n" + "\n")
check_file.write('Number of residues:  ' + str(loaded_mol.num_res()) + "\n")
check_file.write('Number of atoms:  ' + str(loaded_mol.num_atoms()) + "\n" + "\n")
check_file.write('Time:  ' + time.strftime("%a, %d %b %Y %H:%M:%S") + "\n" + "\n")

###CHECK FOR LIGANDS AND UNKNOWN RESIDUES #############################################################################
check_file.write('Ligands and unknown residues:' + "\n")
ligand_dic = loaded_mol.get_ligands_unknown_res()
if ligand_dic == {}:
    print 'There are no ligands or unknown residues in the molecule'
    check_file.write('None' + "\n" + "\n")
else:
    print 'There are ligands or unknown residues.'
    ligand_list = ligand_dic.keys()
    for res_id in ligand_dic.keys():
        print str(res_id) + ' - ' + loaded_mol.res_name_from_id(res_id)
        check_file.write(str(res_id) + ' - ' + loaded_mol.res_name_from_id(res_id) + "\n")
    check_ligand = raw_input('Shall some or all of the above be treated as normal residues, all other residues will'
                             ' be treated as fully atomistic? (y/n) ')
    if check_ligand == ('y' or 'yes'):
        while 1:
            ligand_inp = raw_input('Enter the residue numbers separated by spaces (enter a for all residues): ')
            if ligand_inp == 'a':
                print 'All residues/ligands are treated normally.'
                check_file.write('All are treated normally.' + "\n" + "\n")
                break
            else:
                try:
                    res_selected = [int(s) for s in ligand_inp.split()]
                    new_ligand_list = []
                    string_ligands = ''
                    for i in ligand_list:
                        if i not in res_selected:
                            new_ligand_list.append(i)
                            string_ligands += str(i) + ' '
                    res_atomistic += new_ligand_list
                    print 'The following residues will be treated atomistically: ' + string_ligands
                    check_file.write('The following residues will be treated atomistically: ' + string_ligands
                                     + "\n" + "\n")
                    break
                except ValueError:
                    print 'Please enter integer numbers.'
    else:
        res_atomistic += ligand_list
        print 'The above listed residues/ligands will be treated fully atomistically.'
        check_file.write('All are treated fully atomistically.' + "\n" + "\n")

###SELECTION OF RESIDUES FOR FULLY AND LOCALLY RIGID REGIONS ##########################################################
print
print
print 'Selection of atomistic and locally rigidified residues'
print
check_file.write('Selection of atomistic and locally rigidified residues' + "\n")
res_frozen = remove_selection_from_list(res_dic.keys(), res_atomistic)
x = choose_res(loaded_mol)  #initiates the first round of selections
y = merge_new_and_old_res_lists(res_frozen, res_local, res_atomistic, x[1], x[0])
res_frozen = y[2]
res_local = y[1]
res_atomistic = y[0]
log_string = ''
for strings in x[2]:
    log_string += strings + "\n"
check_file.write(log_string + "\n" + "\n")

while 1:  #allows to choose more residues with different methods
    print
    print 'The following residues are not frozen (atomistic):'
    print res_atomistic
    print
    print 'The following residues are selected to be locally rigidified:'
    print res_local
    print
    print 'The following residues are selected to be fully rigidified:'
    print res_frozen
    print
    check = raw_input('Enter more residues (1) or remove residues from previous selection (2)? Any other key to exit. ')
    try:
        value = int(
            check)  #to enter more residues the choose_res function is called again with the same options as before
        if int(check) == 1:
            x = choose_res(loaded_mol)
            y = merge_new_and_old_res_lists(res_frozen, res_local, res_atomistic, x[1], x[0])
            res_frozen = y[2]
            res_local = y[1]
            res_atomistic = y[0]
            log_string = ''
            for strings in x[2]:
                log_string += strings + "\n"
            check_file.write(log_string + "\n" + "\n")
        elif int(check) == 2:
            x = remove_res(loaded_mol)
            y = removal_selection(res_frozen, res_local, res_atomistic, x[1],
                                  x[0])  #the user enters a selection for removal
            res_frozen = y[2]
            res_local = y[1]
            res_atomistic = y[0]
            check_file.write('Remove residues' + "\n" + 'Atomistic residues removed:' + "\n")
            string = ''
            for i in x[0]:
                string = string + str(i) + '  '
            check_file.write(string + "\n" + "\n")
            check_file.write('Locally rigidified residues removed:' + "\n")
            string = ''
            for i in x[1]:
                string = string + str(i) + '  '
            check_file.write(string + "\n" + "\n")
        else:
            break

    except ValueError:
        break
check_file.write('Selection completed' +  "\n" + "\n")
string_atomistic = ''
for res in res_atomistic:
    string_atomistic += str(res) + ', '
string_local = ''
for res in res_local:
    string_local += str(res) + ', '
check_file.write('Atomistic residues: ' + string_atomistic + "\n" )
check_file.write('Locally rigidified residues: ' + string_local + "\n" )

###GROUPING OF THE RIGID BODIES #######################################################################################
if res_frozen != []:
    for i in res_frozen:
        atoms_used += res_dic[i]
    print
    print 'Grouping of the rigid body'
    print
    print '1-Group all parts of the rigid system as one body'
    print '2-Define groups'
    print
    while 1:
        check = raw_input('Grouping method:  ')
        try:
            value = int(check)
            if int(check) == 1:  #all residues are grouped as one rigid body
                num_rbody += 1
                check_file.write('All frozen residues are rigidified as one body' + "\n" + "\n" + "\n")
                groups = {1: res_frozen}
                break
            elif int(check) == 2:  #residues are frozen in a number of rigid bodies
                check_file.write('The frozen residues are rigidified as more than one body' + "\n" + "\n" + "\n")
                print 'The following residues are frozen:'
                print res_frozen
                print
                res_frozen_left = res_frozen
                num_groups = int(raw_input('How many groups will be defined?  '))  #number of rigid bodies is defined
                counter = 0
                while counter < num_groups:
                    counter += 1
                    num_rbody += 1
                    group = raw_input(
                        'Residues for rigid body: ')  #input allows to enter 'r num1 num2' being the range between num1 and num2
                    group = group.split()
                    i = 0
                    group_res = []
                    while i < len(group):
                        if group[i] == 'r':  #if a range is indicated
                            i += 1
                            start_range = int(group[i])  #finds starting residue
                            i += 1
                            end_range = int(group[i])  #finds final residue
                            for j in range(start_range, end_range):
                                if j in res_frozen_left:  #appends list if residues have not be assigned earlier
                                    group_res.append(j)
                        else:
                            try:
                                if int(group[i]) in res_frozen_left:
                                    group_res.append(int(group[i]))  #appends single residues entered
                                i += 1
                            except ValueError:
                                i += 1
                    res_frozen_left = remove_selection_from_list(res_frozen_left,
                                                                 group_res)  #removes newly assigned residues from rest of residues
                    groups[num_rbody] = group_res  #defines a new entry in the dictionary that saves the groups
                    if res_frozen_left == []:  #automatically terminates when no residues are left to be assigned even if the number of rigid bodies entered is not reached
                        print 'All residues are assigned.'
                        break
                    if counter == num_groups and res_frozen_left != []:  #if all rigid bodies are assigned and there are still residues not assigned it deletes all groups and restarts the process
                        print 'Not all residues were assigned. Please start again.'
                        counter = 0
                        groups = {}
                        res_frozen_left = res_frozen
                break
            else:
                pass
        except ValueError:
            print 'Wrong input'
    print
    print 'The following rigid bodies have been specified:'  #prints all the entered rigid bodies from above
    for i in groups.keys():
        print 'Rigid body: ', i
        print groups[i]
        print

else:
    groups = {}

### LOCALLY RIGIDIFICATION, STANDARD GROUPS ###########################################################################
if res_local != []:
    print 'Grouping of the local rigid bodies'
    print
    local_scheme = []

    proline_check = raw_input('Group proline as rigid body C-O-N-CD-CG-CB-CA(1) or as  N-CD-CG-CB-CA(2)? ')
    if proline_check == '1':
        check_file.write('Grouped proline as rigid body C-O-N-CD-CG-CB-CA' + "\n")
        local_scheme.append(1)
    elif proline_check == '2':
        check_file.write('Grouped proline as rigid body  N-CD-CG-CB-CA' + "\n")
        local_scheme.append(2)
    else:
        local_scheme.append(0)

    arginine_check = raw_input('Group arginine as rigid body NE-HE-CZ-NH1-HH11-HH12-NH2-HH21-HH22 (y/n)? ')
    if arginine_check == 'y':
        check_file.write('Grouped arginine as rigid body NE-HE-CZ-NH1-HH11-HH12-NH2-HH21-HH22' + "\n")
        local_scheme.append(1)
    else:
        local_scheme.append(0)

    histidine_check = raw_input('Group histidine (HIS, HIE, HID, HIP) as rigid body CG-ND1-CE1-NE2-CD2 (y/n)? ')
    if histidine_check == 'y':
        check_file.write('Grouped histidine (HIS, HIE, HID, HIP) as rigid body CG-ND1-CE1-NE2-CD2' + "\n")
        local_scheme.append(1)
        local_scheme.append(1)
        local_scheme.append(1)
        local_scheme.append(1)
    else:
        local_scheme.append(0)
        local_scheme.append(0)
        local_scheme.append(0)
        local_scheme.append(0)

    lysine_check = raw_input('Group lysine as rigid body NZ-HZ1-HZ2-HZ3 (y/n)? ')
    if lysine_check == 'y':
        check_file.write('Grouped lysine as rigid body NZ-HZ1-HZ2-HZ3' + "\n")
        local_scheme.append(1)
    else:
        local_scheme.append(0)

    aspartic_check = raw_input('Group aspartic acid as rigid body CB-CG-OD1-OD2 (y/n)? ')
    if aspartic_check == 'y':
        check_file.write('Grouped aspartic acid as rigid body CB-CG-OD1-OD2' + "\n")
        local_scheme.append(1)
    else:
        local_scheme.append(0)

    asparagine_check = raw_input('Group asparagine as rigid body CB-CG-OD1-ND2-HD21-HD22 (y/n)? ')
    if asparagine_check == 'y':
        check_file.write('Grouped asparagine as rigid body GB-CG-OD1-ND2-HD21-HD22' + "\n")
        local_scheme.append(1)
    else:
        local_scheme.append(0)

    glutamic_check = raw_input('Group glutamic acid as rigid body CG-CD-OE1-OE2 (y/n)? ')
    if glutamic_check == 'y':
        check_file.write('Grouped glutamic acid as rigid body CG-CD-OE1-OE2' + "\n")
        local_scheme.append(1)
    else:
        local_scheme.append(0)

    glutamine_check = raw_input('Group glutamine as rigid body CG-CD-OE1-NE2-HE21-HE22 (y/n)? ')
    if glutamine_check == 'y':
        check_file.write('Grouped glutamine as rigid body CG-CD-OE1-NE2-HE21-HE22' + "\n")
        local_scheme.append(1)
    else:
        local_scheme.append(0)

    phenylalanine_check = raw_input('Group phenylalanine as rigid body CG-CD1-HD1-CE1-HE1-CZ-HZ-CE2-HE2-CD2-HD2 (y/n)? ')
    if phenylalanine_check == 'y':
        check_file.write('Grouped phenylalanine as rigid body CG-CD1-HD1-CE1-HE1-CZ-HZ-CE2-HE2-CD2-HD2' + "\n")
        local_scheme.append(1)
    else:
        local_scheme.append(0)

    tyrosine_check = raw_input('Group tyrosine as rigid body CG-CD1-HD1-CE1-HE1-CZ-OH-CE2-HE2-CD2-HD2 (y/n)? ')
    if tyrosine_check == 'y':
        check_file.write('Grouped tyrosine as rigid body CG-CD1-HD1-CE1-HE1-CZ-OH-CE2-HE2-CD2-HD2' + "\n")
        local_scheme.append(1)
    else:
        local_scheme.append(0)

    tryptophan_check = raw_input(
        'Group tryptophan as either: 1 rigid body:(CG-CD1-HD1-NE1-HE1-CE2-CZ2-HZ2-CH2-HH2-CZ3-HZ3-CE3-HE3-CD2) '
        'or 2 rigid bodies:(CG-CD1-HD1-NE1-HE1)-(CE2-CZ2-HZ2-CH2-HH2-CZ3-HZ3-CE3-HE3-CD2)? ')
    if tryptophan_check == '1':
        check_file.write(
            'Grouped tryptophan as rigid body CG-CD1-HD1-NE1-HE1-CE2-CZ2-HZ2-CH2-HH2-CZ3-HZ3-HZ3-CE3-HE3-CD2' + "\n" + "\n")
        local_scheme.append(1)
    elif tryptophan_check == '2':
        check_file.write(
            'Grouped tryptophan as two rigid bodies CG-CD1-HD1-NE1-HE1 and CE2-CZ2-HZ2-CH2-HH2-CZ3-HZ3-HZ3-CE3-HE3-CD2' + "\n" + "\n")
        local_scheme.append(2)
    else:
        local_scheme.append(0)

    loc = loaded_mol.get_rigid_for_local(res_local, local_scheme, atoms_used)
    groups_local = loc[0]
    atoms_used = loc[1]
    print


    check_peptide_bonds = raw_input('Shall the peptide bonds be treated as rigid bodies (excl. PRO) (y/n)? ')
    if check_peptide_bonds == 'y' or 'yes':
        pairlist = []  #all peptide bonds are part of two residues
        for i in res_local:
            if (i - 1, i) not in pairlist and (i - 1 > 0):
                if loaded_mol.get_residue_name(i) != 'PRO':
                    pairlist.append((i - 1, i))
            if (i, i + 1) not in pairlist and (i + 1 <= loaded_mol.num_res()):
                if loaded_mol.get_residue_name(i + 1) != 'PRO':
                    pairlist.append((i, i + 1))
        for j in pairlist:
            atomlist = []
            atomlist.append(loaded_mol.get_atom_id_from_name('C', j[0]))
            atomlist.append(loaded_mol.get_atom_id_from_name('O', j[0]))
            atomlist.append(loaded_mol.get_atom_id_from_name('N', j[1]))
            atomlist.append(loaded_mol.get_atom_id_from_name('H', j[1]))
            check_double = 0
            for k in atomlist:
                if k in atoms_used:
                    check_double = 1
            if check_double == 0:
                atoms_used += atomlist
                check_file.write(
                    'Peptide bond rigidified for residues ' + str(j[0]) + ' and ' + str(j[1]) + "\n")
                groups_local[len(groups_local) + 1] = [j, (loaded_mol.get_residue_name(j[0]),
                                                           loaded_mol.get_residue_name(j[1]))] + atomlist
    check_file.write("\n")
    check_new_groups = raw_input('Enter additional user-defined rigid bodies (y/n)? ')
    if check_new_groups in ['yes', 'y']:
        print '1 - select a residues by name'
        print '2 - select one or more residues by id'
        while 1:
            while 1:
                local_method = raw_input('Method (1/2): ')
                try:
                    if int(local_method) in [1, 2]:
                        break
                    else:
                        print 'Invalid choice'
                except ValueError:
                    print 'Invalid input'

            if int(local_method) == 1:
                while 1:
                    exit_1 = 0
                    res_name = raw_input('Residue name (three letter code, n to exit): ')
                    if res_name == 'n':
                        exit_1 = 1
                        break
                    try:
                        res_list = loaded_mol.get_all_res_name(res_local, res_name)
                        if res_list != []:
                            check_file.write('Additional rigid bodies for residues: ' + res_name +  "\n")
                            break
                        else:
                            print 'None of the residues for local rigidification is of that type.'
                    except ValueError:
                        print 'Invalid input'
                if exit_1 == 0:
                    all_atoms_res_name = []
                    atoms_used_name = []
                    for i in res_list:
                        all_atoms_res_name += res_dic[i]
                    for j in all_atoms_res_name:
                        if j in atoms_used:
                            atoms_used_name.append(j)
                    if atoms_used_name != []:
                        print 'The following atoms are already in local rigid bodies: '
                        for k in atoms_used_name:
                            print str(k) + ' ' + atom_dic[k][0] + ' ' + atom_dic[k][1] + ' ' + str(atom_dic[k][2])
                    print
                    print 'The residues %s contain the following atoms:' % res_name
                    atoms_in_res = loaded_mol.get_atom_list(res_list[0])
                    for l in atoms_in_res:
                        print atom_dic[l][0]
                    print
                    atom_seq = raw_input('Enter a sequence of atom names separated by spaces: ')
                    for m in res_list:
                        check_res = 0
                        atomlist = res_dic[m]
                        atoms_rbody = []
                        for n in atomlist:
                            if atom_dic[n][0] in atom_seq.split():
                                atoms_rbody.append(n)
                        if len(atoms_rbody) != len(atom_seq.split()):
                            print 'Not all atoms entered are in the residue.'
                            check_file.write('For residue: ' + str(m) + ' - not all entered atoms in residue' + "\n")
                            check_res = 1
                        if len(atoms_rbody) < 3:
                            print 'Rigid bodies need to have at least 3 atoms.'
                            check_file.write('For residue: ' + str(m) + ' - not 3 or more atoms in residue' + "\n")
                            check_res = 1
                        for o in atoms_rbody:
                            if o in atoms_used_name:
                                print 'Atoms cannot belong to two rigid bodies.'
                                check_file.write('For residue: ' + str(m) + ' - atoms cannot belong to two rigid bodies'
                                                 + "\n" + "\n")
                                check_res = 1
                                break
                        if loaded_mol.check_linear(atoms_rbody, tolerance):
                            print 'Linear arrangement of atoms cannot be a rigid body.'
                            check_file.write('For residue: ' + str(m) + ' - linear group entered' + "\n")
                            check_res = 1
                        if check_res == 0:
                            groups_local[len(groups_local) + 1] = [m, loaded_mol.get_residue_name(m)] + atoms_rbody
                            atoms_used += atoms_rbody
                            check_file.write('For residue: ' + str(m) + ' new group: ' + atom_seq + "\n")
                    check_file.write("\n")
            elif int(local_method) == 2:
                while 1:
                    exit_2 = 0
                    res_inp = raw_input('Residue numbers (n to exit): ')
                    if res_inp == 'n':
                        exit_2 = 1
                        break
                    try:
                        res_list_str = res_inp.split()
                        res_list_int = []
                        for res in res_list_str:
                            res_list_int.append(int(res))
                        if res_list_int != []:
                            check_file.write('Additional rigid bodies for residues: ' + res_inp  + "\n")
                            break
                        else:
                            print 'No residues entered'
                    except ValueError:
                        print 'Invalid input'
                if exit_2 == 0:
                    all_atoms_res = []
                    atoms_used_res = []
                    for i in res_list_int:
                        all_atoms_res += res_dic[i]
                    for j in all_atoms_res:
                        if j in atoms_used:
                            atoms_used_res.append(j)
                    atoms_allowed = []
                    for atom in all_atoms_res:
                        if atom not in atoms_used_res:
                            atoms_allowed.append(atom)
                    if atoms_used_res != []:
                        print 'The following atoms are already in local rigid bodies: '
                        for k in atoms_used_res:
                            print str(k) + ' ' + atom_dic[k][0] + ' ' + atom_dic[k][1] + ' ' + str(atom_dic[k][2])
                    print
                    print 'The residues %s contain the following atoms that are not assigned to rigid bodies:' % res_list_str
                    string_atoms = ''
                    for atom in atoms_allowed:
                        string_atoms += str(atom) + ' - ' + atom_dic[atom][0] + ', '
                    print string_atoms
                    print
                    atom_seq = raw_input('Enter a sequence of atom numbers separated by spaces: ')
                    check_res = 0
                    atoms_rbody = []
                    check_rbody = 0
                    for atom_str in atom_seq.split():
                        try:
                            if int(atom_str) in atoms_allowed:
                                atoms_rbody.append(int(atom_str))
                            elif int(atom_str) in atoms_used_res:
                                print 'Atoms cannot belong to two rigid bodies.'
                                check_file.write('Atoms cannot belong to two rigid bodies'  + "\n")
                                check_rbody = 1
                                break
                        except ValueError:
                            pass
                    print len(atoms_rbody)
                    print len(atom_seq.split())
                    if len(atoms_rbody) != len(atom_seq.split()):
                        print 'Not all atoms entered are in the residue.'
                        check_file.write('Not all entered atoms in residue' + "\n")
                        check_rbody = 1
                    if len(atoms_rbody) < 3:
                        print 'Rigid bodies need to have at least 3 atoms.'
                        check_file.write('Not 3 or more atoms in residue' + "\n")
                        check_rbody = 1
                    if loaded_mol.check_linear(atoms_rbody, tolerance):
                        print 'Linear arrangement of atoms cannot be a rigid body.'
                        check_file.write('Linear group entered' + "\n")
                        check_rbody = 1
                    if check_rbody == 0:
                        groups_local[len(groups_local) + 1] = ['', ''] + atoms_rbody
                        atoms_used += atoms_rbody
                        check_file.write('For residues: ' + res_inp + ' new group: ' + atom_seq + "\n")
                    check_file.write("\n")
            check_more = raw_input('Enter more rigid bodies (y/n)? ')
            if check_more not in ['y', 'yes']:
                break

### WRITING OUTPUT FILES ##############################################################################################

print 'Output files will be written now.'

### COORDSINIRIGID ###
coords_out = open('coordsinirigid', 'w')
for i in range(1, len(atom_dic.keys()) + 1):
    atom = str(atom_dic[i][3]) + '   ' + str(atom_dic[i][4]) + '   ' + str(atom_dic[i][5])
    coords_out.write(atom)
    if i < len(atom_dic.keys()):
        coords_out.write("\n")
coords_out.close()

### RBODYCONFIG FILE ###
groups_file = open('rbodyconfig', 'w')
for group in range(1, len(groups.keys()) + 1):
    string = ''
    counter = 0
    for res in groups[group]:
        for atom in res_dic[res]:
            string += str(atom) + "\n"
            counter += 1
    groups_file.write('GROUP ' + str(counter) + "\n" + string)
for group_local in range(1, len(groups_local.keys()) + 1):
    string = ''
    counter = 0
    for atom in groups_local[group_local][2:]:
        string += str(atom) + "\n"
        counter += 1
    groups_file.write('GROUP ' + str(counter) + "\n" + string)
groups_file.close()

### FINISH LOG FILE ###
check_file.close()

### Writing in file parmed.py to exclude interactions ####
parmed_in = open('parmed_in', "w")
for group in range(1, len(groups.keys()) + 1):
    string = 'addExclusions '
    string2 = '@'
    counter = 1
    for res in groups[group]:
	for atom in res_dic[res]:
	    if counter == 1:
                string2 += str(atom)
                counter = 2
            else:
                string2 += ',' + str(atom)
    parmed_in.write(string + string2 + ' ' +string2 +"\n")
for group_local in range(1, len(groups_local.keys()) + 1):
    string = 'addExclusions '
    string2 = '@'
    counter = 1
    for atom in groups_local[group_local][2:]:
        if counter == 1:
            string2 += str(atom)
            counter = 2
        else:
            string2 += ',' + str(atom)
    parmed_in.write(string + string2 + ' ' + string2 + "\n")
parmed_in.write('parmout coords.prmtop.ni' + "\n" + 'go' + "\n")


print 'Selection process completed - change topology file to exclude interactions within the rigid bodies.'

if pymol_check:
    mol_view = PymolView(pdb_inp)
    mol_view.apply_colour_scheme(res_atomistic, res_eical, 2, 3)
    pymol.cmd.save('rigid_%s' % pdb_inp, format='png')
    time.sleep(10)
    mol_view.__del__()
