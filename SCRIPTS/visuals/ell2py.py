#! /usr/bin/env python

"""
Converts a pysites.xyz file and coords.1 file to a pymol script

Scott Olesen
swo24@cam.ac.uk
July 2012

Run in a directory with pysites.xyz and coords.1, or alter the filenames at the bottom of this script. I would add argument parsing using argparse if I were sure that everyone knew what to do with argparse.

To run the script in pymol, use 'pymol -r ell.pml'. By uncommenting the png lines at the bottom of the pymol script, you can tell pymol to produce a png. You can do this in batch by calling 'pymol -c -r ell.pml'.

Colors can be added by adding "color [color_name]" to the end of a line in pysites.xyz. A list of color names native to pymol is at pymolwiki.org/color_values. (I like 'red' and 'marine'.) Different colors can be specified using cmd.set_color("color_name", [r, g, b]) in the script preamble. Color specifications can be left out of pysites.xyz and specified inside pymol.
"""

# add system path to help python find packages while on the cluster
import sys
sys.path.append('/usr/local/lib64/python2.4/site-packages')

import numpy as np
import math

class Molecule:
    """Model PY molecule produced from pysites.xyz"""

    def __init__(self, sites):
        self.sites = sites

    def make_ellipsoids(self, pos, aaxis):
        """Produce a list of ellipsoids corresponding to this molecule at given lab frame coordinates"""
        ells = []
        for site in self.sites:
            ells.append(site.make_ellipsoid(pos, aaxis))

        return ells

class Site:
    """A PY site containing all six semiaxes"""

    def __init__(self, pos, aaxis, semiaxes, color=np.array([1, 0, 0])):
        """
        Create site instance

        pos -- np array of molecule frame xyz position
        aaxis -- np array of molecule frame aaxis coordinates
        semiaxes -- np array of a11, a12, ..., a23
        color -- np array of r, g, b

        """

        self.pos = pos
        self.aaxis = aaxis
        self.semiaxes = semiaxes
        self.color = color

        self.scale = 1.0

        # set up molecule frame semiaxis vectors from repulsive ellipsoid
        self.rotmat = aaxis_to_mat(aaxis)
        self.e = []
        self.e.append(np.dot(self.rotmat, np.array([semiaxes[0], 0.0, 0.0])))
        self.e.append(np.dot(self.rotmat, np.array([0.0, semiaxes[1], 0.0])))
        self.e.append(np.dot(self.rotmat, np.array([0.0, 0.0, semiaxes[2]])))

    def make_ellipsoid(self, pos, aaxis):
        rotmat = aaxis_to_mat(aaxis)
        ell_pos = pos + np.dot(rotmat, self.pos)
        ell_scale = self.scale
        ell_e = []
        for ei in self.e:
            ell_e.append(np.dot(rotmat, ei))
        ell_color = self.color

        return Ellipsoid(ell_pos, ell_scale, ell_e, ell_color)

class Ellipsoid:
    """An ellipsoid, used for visualizing PY sites"""

    def __init__(self, pos, scale, e, color):
        self.pos = pos
        self.scale = scale
        self.e = e
        self.color = color

    def cgo_string(self, name):
        """Produce CGO string"""

        # cgo.ELLIPSOID, x0, y0, z0, scale, a1, a2, a3, b1, b2, b3, c1, c2, c3,
        s = 'cmd.load_cgo([cgo.ELLIPSOID,'
        for x in self.pos:
            s = s + str(x) + ','
        s = s + str(self.scale) + ','
        for ei in self.e:
            for x in ei:
                s = s + str(x) + ','
        s = s[:-1] + '], "{0}")\n'.format(name)

        # set color
        if self.color is not None:
            s = s + 'cmd.color("{0}", "{1}")\n'.format(self.color, name)

        return s


def aaxis_to_mat(p):
    """Converts an angle-axis rotation into a rotation matrix"""
    # cf wikipedia page for Rodrigues's rotation formula
    theta = np.linalg.norm(p)
    if theta == 0:
        return np.identity(3)
    else:
        k = p / theta
        kx = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
        return np.identity(3) * math.cos(theta) + kx * math.sin(theta) + (1 - math.cos(theta)) * np.outer(k, k)

def xyz_to_molecule(lines):
    """
    Produce an array of Molecule instances from pysites.xyz

    lines -- pysites.xyz split into lines

    """

    s = [line.split() for line in lines]

    # check that the pysites file was properly written
    nsites = int(s[0][0])
    assert(nsites + 2 == len(lines))

    sites = []
    for i in range(2, nsites + 2):
        pos = np.array([float(x) for x in s[i][1:4]])
        semiaxes = np.array([float(x) for x in s[i][7:13]])
        aaxis = np.array([float(x) for x in s[i][14:17]])

        # check for color specification in pysites.xyz
        if len(s[i]) == 19 and s[i][17] == 'color':
            color = s[i][18]
        else:
            color = None

        sites.append(Site(pos, aaxis, semiaxes, color))

    return Molecule(sites)

def pysites_to_molecule(filename):
    """Convert pysites.xyz file to a model molecule"""
    with open(filename, 'r') as f:
        lines = f.readlines()

    return xyz_to_molecule(lines)

def coords_to_ellipsoids(filename, model):
    """Convert coords.1 file to ellipsoid list using model molecule"""

    # read in the coords file
    with open(filename, 'r') as f:
        coords = [[float(x) for x in line.split()] for line in f.readlines()]

    # rearrange the coords information
    # pos1, pos2, ..., aax1, aax2, ... -> (pos1, aax1), (pos2, aax2), ...
    coords = zip(coords[:len(coords)/2], coords[len(coords)/2:])

    # produce ellipsoids from the coords
    ells = []
    for (pos, aaxis) in coords:
        ells = ells + model.make_ellipsoids(pos, aaxis)
    
    return ells

def ells_to_script(filename, ells):
    """Convert the ellipsoid list into a pymol script"""

    # produce the pymol string
    s = \
"""import cgo
cmd.reset()
cmd.bg_color("white")
cmd.set("depth_cue", 0)
cmd.set("ellipsoid_quality", 3)
cmd.set("gamma", 1.2)
cmd.set("ambient", 0.2)
cmd.set("spec_reflect", 0.9)
"""

    for i in range(len(ells)):
        s = s + ells[i].cgo_string('ell{0}'.format(i))

    # add ray tracing and rendering lines
    s = s + \
"""#cmd.ray(100,100)
#cmd.png('test', dpi=300)"""

    # save string to file
    with open(filename, 'w') as f:
        f.write(s)
    

if __name__ == '__main__':
    # read the pysites file and produce a model molecule
    model = pysites_to_molecule('pysites.xyz')

    # read coords and produce ellipsoids
    ells = coords_to_ellipsoids('coords.1', model)

    # produce pymol script
    ells_to_script('ell.py', ells)
