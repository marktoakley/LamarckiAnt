#!/bin/bash
# COMPONENT SCRIPT FOR ITERATIVE RESP CHARGE GENERATION (WRITE RESP CHARGES TO NEW PRMTOP)
# Original script by Kyle Sutherland-Cash, modified by Chris Whittleston (csw34)

#########################################################################################
# You shouldn't need to modify anything in this script as it is called by run_gamess.sh #
#########################################################################################

# Writes an Antechamber format file given a PDB ($1) and net charge ($2)
function write_ac_file {
    # Give it a PDB to work out file names from
    pdb=$1
    file_prefix=`echo $1 | sed 's/\.pdb//g'`

    # Assign other arguments 
    gamess_log=$2
    charge=$3

    # Write a charge file containing zeroes for every line
    sed 's/.*/0/g' ${pdb} > ./tmp.cf

    # Write everything preceding the coordinates (columns 31-54 inclusive) of the PDB to a file
    # and everything after the coordinates too
    awk '{print substr($0, 1, 30);}' ${pdb} > ./tmp.start
    awk '{print substr($0, 55);}' ${pdb} > ./tmp.end

    # Write the coordinates in Bohr to a coordinates file from the gamess_log file
    sed -n '/COORDINATES (BOHR)/,/ATOMIC BASIS SET/ p' ${gamess_log} | \
    awk 'BEGIN {bohr = 0.5291772108; }
         /^ [HCNOS]/ {printf "%8.3f%8.3f%8.3f\n", $3 * bohr, $4 * bohr, $5 * bohr;}' > ./tmp.coords

    # Paste them together to make ./tmp.pdb
    paste ./tmp.start ./tmp.coords ./tmp.end > ./tmp.pdb
    rm ./tmp.start ./tmp.coords ./tmp.end

    # Use Antechamber to write to ${file_prefix}.tmp.ac
    antechamber -i ./tmp.pdb \
                -fi pdb \
                -o ./tmp.ac \
                -fo ac \
                -c rc \
                -nc ${charge} \
                -cf ./tmp.cf \
                -pf y > /dev/null
    # Remove temporary files
    rm ./tmp.pdb ./tmp.cf
}

# Writes an ESP file given a GAMESS log file ($1)
function write_esp_file {
    # Assign arguments
    gamess_log=$1

    # Write the number of atoms and points for fitting and zero to the top of the ESP file
    awk '/TOTAL NUMBER OF ATOMS/ {natoms = $6;} 
         /NUMBER OF POINTS SELECTED FOR FITTING/ {nesp = $8;}
         END {printf "%5d%5d%5d\n", natoms, nesp, 0;}' ${gamess_log}

    # Write the coordinates in Bohr to the ESP file from the gamess_log file
    sed -n '/COORDINATES (BOHR)/,/ATOMIC BASIS SET/ p' ${gamess_log} | \
    awk '/^ [HCNOS]/ {printf "%33.7E%16.7E%16.7E\n", $3, $4, $5;}' 
     

    # Echo the log file to sed and select the lines referring to electrostatic potential
    sed -n '1,/NUMBER OF POINTS SELECTED FOR FITTING/ !p' ${gamess_log} | \
    sed -n '/END OF PROPERTY EVALUATION/,$ !p' | \
    sed 's/^.\{5\}//g' | \
    awk '{printf "%17.7E%16.7E%16.7E%16.7E\n", $6, $1, $2, $3;}'
}

# Function to generate a qout file from tmp.ac ($1) and tmp.esp ($2)
function generate_qout {
    # Assign arguments
    ac_file=$1
    esp_file=$2

    # Clear temporary files in case they already exist
    rm ./tmp.respout[12] ./tmp.respin[12] ./tmp.qout[12] ./esout ./punch 2> /dev/null

    # Run respgen to generate respin files for resp charge calculation
    respgen -i ${ac_file} \
            -o ./tmp.respin1 \
            -f resp1 2>&1 1> /dev/null
    respgen -i ${ac_file} \
            -o ./tmp.respin2 \
            -f resp2 2>&1 1> /dev/null

    # Run resp to calculate the charges and remove the temporary files
    resp -i ./tmp.respin1 \
         -o ./tmp.respout1 \
         -e ${esp_file} \
         -t ./tmp.qout1 2>&1 1> /dev/null
    rm ./punch ./esout 2>&1 1> /dev/null

    # Run resp a second time and remove the temporary files again
    resp -i ./tmp.respin2 \
         -o ./tmp.respout2 \
         -e ${esp_file} \
         -q ./tmp.qout1 \
         -t ./tmp.qout2 2>&1 1> /dev/null

    # Print them out in the same format as the topology file (conversion factor of 18.2223 for Amber prmtop)
    cat ./tmp.qout2 | \
    awk 'BEGIN { j = 0
                 printf "%-80s\n", "%FORMAT(5E16.8)"
               }
         { 
            for (i = 1; i <= NF; i++) {
                j++
                printf "%16.8E", $i * 18.2223
                if (j % 5 == 0) {
                    printf "\n"
                }
            }
         }
         END { printf "\n" }'

    # Remove temporary files
    rm ./tmp.respin[12] ./tmp.respout[12] ./tmp.qout[12] ${ac_file} ${esp_file} ./esout ./punch 2>&1 1> /dev/null
}

# Assign command line arguments
pdb=$1
gamess_log=$2
charge=$3
old_prmtop=$4
new_prmtop=$5

# Call functions
write_ac_file ${pdb} ${gamess_log} ${charge} > ./tmp.ac
write_esp_file ${gamess_log} > ./tmp.esp
sed -n '1,/%FLAG CHARGE/ p' ${old_prmtop} > ${new_prmtop}
generate_qout ./tmp.ac ./tmp.esp >> ${new_prmtop}
sed -n '/%FLAG MASS/,$ p' ${old_prmtop} >> ${new_prmtop}
