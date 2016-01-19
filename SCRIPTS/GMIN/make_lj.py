import sys
import math
import random

def rand_about_zero(radius):
    return random.random(-radius, radius)

number = int(sys.argv[1])
density = 0.74
volume = number * density
radius = math.pow((volume * 3.0 / (4.0 * math.pi)), 1.0/3.0)

for i in xrange(0, number):
    x, y, z = (random.uniform(-radius, radius) for j in xrange(0, 3))
    print "{:20.10f}{:20.10f}{:20.10f}".format(x, y, z)
