#!/bin/bash

numb=`wc min.data | awk '{print$1}'`

echo $numb 
rm list_all
for  ((i=1; i<=$numb; i+=1))
  do
   echo $i out of $numb
   cat ORDER_PARAM/pathsample.out.$i | grep AMHQ | grep -v SETUP >> list_all
done #  for  ((i=1; i<=12; i+=1))

echo "#!/usr/bin/python " > histo.py
echo " " >> histo.py
echo "import numpy " >> histo.py
echo "import sys " >> histo.py
echo "import getopt " >> histo.py
echo " " >> histo.py
echo "x=[] " >> histo.py
echo "y=[] " >> histo.py
echo " " >> histo.py
echo "filename = sys.argv[1] " >> histo.py
echo " " >> histo.py
echo "f=open(filename,\"r\") " >> histo.py
echo " " >> histo.py
echo "for line in f: " >> histo.py
echo "        tokens=line.split() " >> histo.py
echo "        x.append(float(tokens[1])) " >> histo.py
echo "        y.append(int(tokens[2])) " >> histo.py
echo " " >> histo.py
echo "bins = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0] " >> histo.py
echo " " >> histo.py
echo "total=[] " >> histo.py
echo " " >> histo.py
echo "for a, b in enumerate(bins): " >> histo.py
#echo "#    print "VALUE INDEX %f, %d,  %d " % (b,a,a+1) " >> histo.py    
echo "    total.append([]) " >> histo.py
echo "    for i in range(len(x)): " >> histo.py
echo "       if x[i] >= b and x[i] < bins[a+1] : " >> histo.py
echo "               total[a].append(y[i]) " >> histo.py
echo " " >> histo.py
echo "for a, b in enumerate(bins): " >> histo.py
echo "      filefile=\"data.\"+str(a) " >> histo.py
echo "      numpy.savetxt(filefile, total[a], fmt='%d') " >> histo.py

chmod +x histo.py

echo "DELTA 10.0" > dinfo
echo "FIRST 90.0" >> dinfo
echo "LEVELS 40"  >> dinfo
echo "COMMENT IDMIN 1"  >> dinfo
echo "COMMENT IDMIN 2"  >> dinfo
echo "MINIMA min.data"  >> dinfo
echo "TS ts.data"  >> dinfo
echo "CENTREGMIN"  >> dinfo
echo "TRMIN 5 $numb  data.0 data.1 data.2 data.3 data.4" >> dinfo
echo "COMMENT LOWEST 20000" >> dinfo
echo "COMMENT NCONNMIN 0"  >> dinfo
echo "COMMENT NOBARRIERS"  >> dinfo
echo "COMMENT CONNECTMIN 2498"  >> dinfo

