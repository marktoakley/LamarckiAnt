#!/usr/bin/python
# COMPONENT SCRIPT FOR ITERATIVE RESP CHARGE GENERATION (WRITE GAMESS INPUT FILE)
# Original script by Kyle Sutherland-Cash, modified by Chris Whittleston (csw34)

# CHECK: Make sure ~csw34/gamess/rungms exists on the machine you are running this on!
#        If not - modify the line near the bottom of this file!

#########################################################################################
# You shouldn't need to modify anything in this script as it is called by run_gamess.sh #
#########################################################################################

import re
import sys

def write_output(input_pdb, gamess_file, output_pdb, output_job, charge):
    input = open(input_pdb, 'r')
    output = open(gamess_file, 'w')

    nuc_charges = {'H': '1.0', 'C': '6.0', 'N': '7.0', 'O': '8.0', 'S': '16.0'}
    # Write GAMESS US input file, N.B. GAMESS requires ' $' to be the first
    # two columns of any namelist group.  Otherwise, the formatting is fairly
    # free.  GAMESS will ignore anything outside a $KEYWORD...$END.
    output.write(' $CONTRL  ICHARG=' + charge + ' MPLEVL=0\n')
    output.write('          RUNTYP=ENERGY SCFTYP=RHF MULT=1\n')
    output.write('          MAXIT=200 COORD=CART              $END\n')
    output.write(' $SCF     CONV=1.0E-08 DIRSCF=.T.           $END\n')
    output.write(' $SYSTEM  TIMLIM=50000 MWORDS=120 MEMDDI=0  $END\n')
    output.write(' $BASIS   GBASIS=N31 NGAUSS=6 NDFUNC=1\n')
    output.write('          NPFUNC=0 DIFFSP=.F.               $END\n')
    output.write(' $GUESS   GUESS=HUCKEL                      $END\n')
    output.write(' $ELPOT   IEPOT=1 WHERE=PDC OUTPUT=BOTH     $END\n')
    output.write(' $PDC     PTSEL=CONNOLLY CONSTR=NONE        $END\n')
    output.write(' $DATA\n')
    output.write(' ' + re.split('\/', output_pdb)[-1] + '\n')
    output.write(' C1\n')
    for line in input:
        if 'ATOM' in line:
            x = float(line[30:37])
            y = float(line[38:45])
            z = float(line[46:53])

            atom_name = line[12:15]

            while atom_name[0] not in nuc_charges:
                atom_name = atom_name[1:]

            output.write(' ' + atom_name[0] + '   ')
            output.write(nuc_charges[atom_name[0]] + '   ')
            output.write('%8.3f' % x)
            output.write('%8.3f' % y)
            output.write('%8.3f' % z)
            output.write('\n')
    output.write(' $END')

    input.close()
    output.close()

    gamess_qsub = open(output_job, 'w')

    prefix = re.split('\.', re.split('\/', input_pdb)[-1])[-2]

    gamess_qsub.write('#PBS -q s8\n\n')
    gamess_qsub.write('cd $PBS_O_WORKDIR\n\n')
    gamess_qsub.write('~csw34/gamess/rungms ' + prefix + ' 8 > ' + prefix + '.log')

    gamess_qsub.close()

write_output(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
