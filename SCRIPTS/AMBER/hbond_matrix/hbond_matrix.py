#!/usr/bin/python

def make_hbond_input(hbond_input_filename, hbond_output_filename,
                     donor_acceptor_filename, restart_filename,
                     distance_cutoff, angle_cutoff):
    """
    Make the hbond input file with specified cutoffs and input/output files.
    """
    with open(hbond_input_filename, 'w') as hbond_input:
        # trajin command to read restart file
        hbond_input.write(' '.join(['trajin', restart_filename, '\n']))
        # contents of donor/acceptor file to specify hbonding groups
        with open(donor_acceptor_filename, 'r') as donor_acceptor_file:
            for line in donor_acceptor_file:
                hbond_input.write(line)
        # trajout command to write output file
        hbond_input.write(' '.join(['hbond', 'out', hbond_output_filename,
                                    'distance', str(distance_cutoff),
                                    'angle', str(angle_cutoff), '\n']))

def run_ptraj(hbond_input_filename, topology_filename):
    """
    Run ptraj as a subprocess, using the hbond input file as STDIN and piping
    STDOUT to /dev/null.
    """
    import os
    import os.path
    import subprocess
    with open(hbond_input_filename, 'r') as hbond_input:
        with open(os.devnull, 'w') as devnull:
            # Runs ptraj
            ptraj = subprocess.Popen([os.path.join(os.environ['AMBERHOME'], 'bin', 'ptraj'), topology_filename], stdin=hbond_input, stdout=devnull)
            ptraj.wait()

def get_residues(residues_filename, topology_filename):
    """
    Get the residues from the residue file, if specified. If the file is
    not specified, use the list of all residues from the topology file.
    """
    residues = []
    # If residue file is specified, use those residues.
    if residues_filename:
        with open(residues_filename, 'r') as residue_input:
            for line in residue_input:
                residues.append(int(line.strip()))
    # Otherwise, use all residues from the topology file.
    else:
        with open(topology_filename, 'r') as prmtop:
            target_line_num = None
            for i, line in enumerate(prmtop):
                if 'POINTERS' in line:
                    target_line_num = i + 3
                if i == target_line_num:
                    num_residues = int(line.split()[1])
                    residues = range(1, num_residues + 1)
    return residues

def create_hbond_dict(hbond_output_filename):
    """
    Creates a dictionary containing the hydrogen bonds found by ptraj. This has
    keys which are tuples of the form:

    (donor_res, 'bb'|'sc', acc_res, 'bb'|'sc')
    
    The values are the number of bonds.
    
    e.g.
        hbond_dict[(12, 'bb', 15, 'sc')] = 2 indicates that residue 12 is a 
        backbone donor to the sidechain of residue 15.
    """
    hbond_dict = {}
    with open(hbond_output_filename, 'r') as hbond_output:
        for line in hbond_output:
            words = line.split()
            try:
                if words[0] == '|':
                    donor_res, donor_atom = map(str.strip, words[2].strip(':').split('@'))
                    acceptor_res, acceptor_atom = map(str.strip, words[7].strip(':').split('@'))
                    # Make sure that the residue numbers are integers.
                    donor_res, acceptor_res = map(int, (donor_res, acceptor_res))
                    if donor_atom == 'O' or donor_atom == 'OXT':
                        donor_type = 'bb'
                    else:
                        donor_type = 'sc'
                    if acceptor_atom == 'N':
                        acceptor_type = 'bb'
                    else:
                        acceptor_type = 'sc'
                    try:
                        hbond_dict[(donor_res, donor_type, acceptor_res, acceptor_type)] += 1
                    except KeyError:
                        hbond_dict[(donor_res, donor_type, acceptor_res, acceptor_type)] = 1
            except IndexError:
                continue
    return hbond_dict

def print_matrix(hbond_dict, residues):
    """
    Prints the hbond dictionary as a matrix with two rows per donor residue
    (one backbone, one sidechain) and two columns per acceptor residue.
    """
    for res1 in residues:
        for res1_type in ['bb', 'sc']:
            for res2 in residues:
                for res2_type in ['bb', 'sc']:
                    try:
                        print hbond_dict[(res1, res1_type, res2, res2_type)],
                    except KeyError:
                        print "0",
            print ""

if __name__ == '__main__':
    import sys

    # Assign the command-line arguments to specify file names
    topology_filename, restart_filename, donor_acceptor_filename = sys.argv[1:4]
    # If specified, get a list of residues we're interested in
    if len(sys.argv) >= 5:
        residues_filename = sys.argv[4]
    else:
        residues_filename = None
    # If specified, set a custom distance and angle cutoff
    if len(sys.argv) == 7:
        distance_cutoff, angle_cutoff = map(float, sys.argv[5:7])
    # Otherwise use defaults of 3.0A and 120 degrees.
    else:
        distance_cutoff, angle_cutoff = 3.0, 120.0

    # Set default names for the hbond files
    hbond_input_filename, hbond_output_filename = 'hbond.in', 'hbond.out'

    # Make the input file
    make_hbond_input(hbond_input_filename, hbond_output_filename,
                     donor_acceptor_filename, restart_filename,
                     distance_cutoff, angle_cutoff)

    # Run ptraj
    run_ptraj(hbond_input_filename, topology_filename)

    # Get the list of important residues
    residues = get_residues(residues_filename, topology_filename)

    # Create the dictionary corresponding to the hbond matrix
    hbond_dict = create_hbond_dict(hbond_output_filename)

    # Print the matrix
    print_matrix(hbond_dict, residues)
