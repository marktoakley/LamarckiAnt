#!/bin/bash

#   SETUP script for PATHSAMPLE
#   move params proteins and match directories
#   move gamma.dat pro.list input_amh
#   create odata.connect odata.start odata.finish
#   create endpoint directories for GMIN

PROTNAME=1uzcz

echo $PROTNAME

cp ~/amh/md_input/pro.list/pro.list.$PROTNAME pro.list
cp ~/amh/md_input/input/$PROTNAME/anneal/wt/amc/input1 input_amh
cp ~/amh/md_input/gammas/gammas.mike gamma.dat

mkdir proteins
for j in `cat pro.list | awk '{print$1}'`; do cp ~/amh/proteins/$j proteins ; done

mkdir match
cp -r ~/amh/match/$PROTNAME match

for j in `cat pro.list | awk '{print$1}'`; do cp ~/amh/match/$PROTNAME/$j match/$PROTNAME ; done

cp -r ~/amh/params .

echo 'REOPTIMISEENDPOINTS'  > odata.connect
echo 'NEBK 100.0' >> odata.connect
echo 'MAXSTEP 0.1' >> odata.connect
echo 'MAXMAX 0.2' >> odata.connect
echo 'TRAD 0.2' >> odata.connect
echo 'DIJKSTRA    2' >> odata.connect
echo 'NEWCONNECT 1 1 10.0 25.0 60 3.0 0.001' >> odata.connect
echo 'BFGSTS 400 5 100 0.001' >> odata.connect
echo 'BFGSMIN   1.0D-7' >> odata.connect
echo 'NOHESS' >> odata.connect
echo 'ENDHESS' >> odata.connect
echo 'ENDNUMHESS' >> odata.connect
echo 'DUMPALLPATHS' >> odata.connect
echo 'MAXBFGS  0.5' >> odata.connect
echo 'UPDATES 200' >> odata.connect
echo 'PUSHOFF 0.01' >> odata.connect
echo 'STEPS     1000' >> odata.connect
echo 'BFGSSTEPS 30000' >> odata.connect
echo 'AMH' >> odata.connect


echo 'DUMPDATA' > odata.start
echo 'ENDNUMHESS' >> odata.start
echo 'REOPTIMISEENDPOINTS' >> odata.start
echo 'MAXBFGS 0.4' >> odata.start
echo 'BFGSMIN 1.0D-7' >> odata.start
echo 'UPDATES 200' >> odata.start
echo 'BFGSSTEPS 100000' >> odata.start
echo 'AMH' >> odata.start

cp odata.start odata.finish

# setup for endpoint calculations with GMIN

mkdir A
cp -r params match proteins A
cp input_amh gamma.dat pro.list A
echo 'SLOPPYCONV  1.0D-7'   > data
echo 'TIGHTCONV   1.0D-7'   >> data
echo 'UPDATES     500'   >> data
echo 'AMH'   >> data
echo 'NINT_AMH    5'   >> data
echo 'EDIFF 0.001'   >> data
echo 'STEPS  2 1.0'   >> data
echo 'MAXERISE 0.001'   >> data
echo 'STEP  0.25 0.0'   >> data
echo 'MAXIT 50000 500'   >> data
echo 'MAXBFGS   0.10'   >> data
echo 'TEMPERATURE 0.1'   >> data
echo 'RANSEED 1'   >> data
echo 'RADIUS 1000000.0'   >> data

cp data A

mkdir B
cp -r params match proteins B
cp input_amh gamma.dat pro.list B
cp data B

echo '#!/bin/bash'   > qsub.test
echo '#PBS -q l16'   >> qsub.test
echo '#PBS -j oe'   >> qsub.test
echo '#PBS -N AMH_v'   >> qsub.test
echo '##PBS -l walltime=24:00:00'   >> qsub.test
echo '# Needed for clust, not mek-quake'   >> qsub.test
echo '#PBS -W x=NACCESSPOLICY:SINGLEJOB'   >> qsub.test
echo '#PBS -l nodes=16'   >> qsub.test
echo ''   >> qsub.test
echo 'cd $PBS_O_WORKDIR'   >> qsub.test
echo 'cat $PBS_NODEFILE >& output'   >> qsub.test
echo 'wc output > nodes.info'   >> qsub.test
echo 'cat $PBS_NODEFILE >> nodes.info'   >> qsub.test
echo 'echo $USER >> nodes.info'   >> qsub.test
echo 'echo $PBS_O_WORKDIR >> nodes.info'   >> qsub.test
echo ''   >> qsub.test
echo '/home/mp466/svn/PATHSAMPLE/bin/pathsample.2.1 >& output'   >> qsub.test
echo ''   >> qsub.test
echo 'echo'   >> qsub.test
echo 'qstat -f ${PBS_JOBID}@clust'   >> qsub.test
echo 'echo'   >> qsub.test
echo ''   >> qsub.test

echo 'DEBUG'  > pathdata
echo 'AMH'    >>  pathdata
echo 'COPYFILES gamma.dat params match proteins pro.list input_amh start finish'    >>  pathdata
echo 'COPYOPTIM'    >>  pathdata
echo 'NATOMS         189'    >>  pathdata
echo 'SYSTEM         AMH'    >>  pathdata
echo 'TEMPERATURE    1.0'    >>  pathdata
echo 'CONNECTIONS    1'    >>  pathdata
echo 'SEED           1'    >>  pathdata
echo 'PERTURB        0.40'    >>  pathdata
echo 'ETOL           1.0D-7'    >>  pathdata
echo 'GEOMDIFFTOL    0.1'    >>  pathdata
echo 'ITOL           1.0D0'    >>  pathdata
echo 'DIRECTION      AB'    >>  pathdata
echo 'EXEC            /home/mp466/svn/OPTIM/bin/nag/AMHOPTIM.4.0'    >>  pathdata
echo 'TRIPLES'    >>  pathdata
echo 'ADDTRIPLES'    >>  pathdata
echo 'STARTTRIPLES'    >>  pathdata
echo 'COMMENT DIJINITSTART  2'    >>  pathdata
echo 'COMMENT CYCLES        1'    >>  pathdata
echo ''    >>  pathdata
echo 'DIJINITCONTFLY       2'    >>  pathdata
echo 'CYCLES               400'    >>  pathdata
echo ''    >>  pathdata
echo 'COMMENT SHORTCUT     50    BARRIER'    >>  pathdata
echo 'COMMENT CYCLES       25'    >>  pathdata
echo 'COMMENT FREEPAIRS 10.0   20.0'    >>  pathdata
echo 'COMMENT PLANCK  9.5D-14'    >>  pathdata
echo ''    >>  pathdata
echo 'JOBSPERNODE 1'    >>  pathdata
echo 'PAIRLIST 1'    >>  pathdata

