#! /usr/bin/env python

"""
Converts rigid body coordinates to pymol script producing cylinders

Scott Olesen (swo24@cam.ac.uk) August 2012

This script imports argparse. If argparse is not installed in python, you can put argparse.py in the same directory as this script. (hg.python.org/cpython/file/default/Lib/argparse.py)

Run in a directory with chiro.xyz output from GMIN. The script reads the length of the rods from the xyz file. If the simulation was in 2D with simple LJ sites, then a non-physical length must be specified for the visualization. The default can be changed using the -z flag on the command line.

To run the script in pymol, use 'pymol -r chir.py'. Altering and uncommenting the ray and png commands will be helpful when producing multiple images for publication figures, but it might have to be combined with move commands.

There are a lot of ways to play with the visualization. Here I modify the ray trace mode, shadows, ambient lighting, specular reflectivity, gamma, and the depth cue.
"""

# cluster package locations
import sys
sys.path.append('/usr/local/lib64/python2.4/site-packages')

import numpy as np
import angleaxis as aa

class Cylinder:
    def __init__(self, pos, vec, radius, length, color):
        self.pos = pos
        self.vec = vec
        self.radius = radius
        self.length = length
        self.color = color

        # set up end point
        self.end = self.pos + self.length * self.vec

    def cgo_string(self):
        """Produce cgo string"""

        s = 'cgo.CYLINDER,'
        s = add_csv(s, self.pos)
        s = add_csv(s, self.end)
        s = add_csv(s, [self.radius])
        s = add_csv(s, self.color)
        s = add_csv(s, self.color)
        s = s[:-2]  # chomp comma and trailing space
        return s

    def cgo_sphere_string(self):
        s = 'cgo.COLOR,'
        s = add_csv(s, self.color)
        s = s + 'cgo.SPHERE,'
        s = add_csv(s, self.end)
        s = add_csv(s, [self.radius])
        s = s[:-2]
        return s


def add_csv(string, list):
    out = string
    for x in list:
        out = out + str(x) + ', '
    return out

def cyl_to_line(cyl, name):
    return 'cmd.load_cgo([' + cyl.cgo_string + '], "' + name + '"'

def chiro_to_cylinders(file, radius, zerolength=0.0):
    #north_color = [1.0, 1.0, 1.0]
    north_color = np.array([1.0, 1.0, 1.0])
    #south_color = [0.5, 0.5, 0.5]
    south_color = 0.3 * north_color

    with open(file, 'r') as f:
        ls = f.readlines()

    s = [l.split() for l in ls]

    # check that the file was properly written
    nsites = int(s[0][0])

    length = float(s[1][12])
    if length == 0.0:
        length = zerolength

    cyls = []
    for i in range(2, nsites + 2):
        pos = np.array([float(x) for x in s[i][1:4]])
        vec = np.array([float(x) for x in s[i][5:8]])
        cyls.append(Cylinder(pos, vec, radius, length, north_color))
        cyls.append(Cylinder(pos, vec, radius, -length, south_color))

    return cyls

def cyls_to_script(filename, cyls):
    # produce the pymol string
    s = \
"""import cgo
cmd.reset()
cmd.bg_color("white")
cmd.set("depth_cue", 1)
cmd.set("gamma", 1.2)
cmd.set("ambient", 0.3)
cmd.set("spec_reflect", 0.8)
cmd.set("ray_trace_mode", 1)
cmd.set("ray_shadows", 0)
"""

    for i in range(0, len(cyls), 2):
        s = s + 'cmd.load_cgo([{0}], "north_{1}")\n'.format(cyls[i].cgo_string(), str(i))
        s = s + 'cmd.load_cgo([{0}], "north_sphere_{1}")\n'.format(cyls[i].cgo_sphere_string(), str(i))
        s = s + 'cmd.load_cgo([{0}], "south_{1}")\n'.format(cyls[i + 1].cgo_string(), str(i + 1))
        s = s + 'cmd.load_cgo([{0}], "south_sphere_{1}")\n'.format(cyls[i + 1].cgo_sphere_string(), str(i + 1))

    # add ray tracing and rendering lines
    s = s + \
"""cmd.reset()
#cmd.ray(100,100)
#cmd.png('test', dpi=300)"""

    # save string to file
    with open(filename, 'w') as f:
        f.write(s)
    

if __name__ == '__main__':
    import argparse

    # Set up argument parsing
    parser = argparse.ArgumentParser(description='Convert chiro.xyz to pymol script')
    parser.add_argument('-f', '--filename', type=str, default='chiro.xyz')
    parser.add_argument('-o', '--output', type=str, default='chiro.py')
    parser.add_argument('-r', '--radius', type=float, default=0.1)
    parser.add_argument('-z', '--zerolength', type=float, default=1.0, help='Length of representational rod for LJ sites')
    rawargs = vars(parser.parse_args()).items()
    args = dict(filter(lambda item: item[1] is not None, rawargs))

    # read the chiro.xyz file and produce a list of ellipsoids
    cyls = chiro_to_cylinders(args['filename'], args['radius'], args['zerolength'])

    # produce pymol script
    cyls_to_script(args['output'], cyls)
