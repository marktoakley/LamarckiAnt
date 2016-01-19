#!/usr/bin/python

import glob
import collections as coll
sub_defined_in = coll.defaultdict(list)
sub_called_in = coll.defaultdict(list)
subs_defined = coll.defaultdict(list)
subs_called = coll.defaultdict(list)

for source_filename in glob.glob('./*.[fF]*'):
   with open(source_filename, 'r') as source_file:
      for line in source_file:
         if len(line.lower().lstrip()) == 0:
            continue
         if line.lower().lstrip()[0] == '!' or line.lower().lstrip()[0] == 'c':
            continue
         if line.lower().lstrip()[0:3] == 'end':
            continue
         words = line.lower().split()
         if 'subroutine' in words[0]:
            try:
               sub_name = words[1].split('(')[0]
            except IndexError:
               print "Oh no:", words
               raise Exception
            sub_defined_in[sub_name] = source_filename
            subs_defined[source_filename] += [sub_name]
         if 'call' in words:
            if 'call' != words[0] and 'if' != words[0]:
               continue
            for i, word in enumerate(words):
               if word == 'call':
                  try:
                     sub_name = words[i+1].split('(')[0]
                  except IndexError:
                     print "Oh no:", words
            sub_called_in[sub_name] = source_filename
            subs_called[source_filename] += [sub_name]

for source_file in subs_called:
   print "================================"
   print source_file
   print "********************************"
   for sub_called in subs_called[source_file]:
      print "calls:", sub_called
      print "from:", sub_defined_in[sub_called]
   print "================================"
