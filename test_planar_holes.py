"""
    description: a python 2.7 program to sample n points u.i.r. 
             from an integer grid or a disk.
             the program approximates the leading constant of 
             the quadratic term of the number of k-holes 
             and outputs the mean value and the standard deviation.
    (c) 2020 Manfred Scheucher <scheucher@math.tu-berlin.de>
"""

from itertools import *
from sys import * 
from random import *
from math import *
import PyDCG # download from http://rfabila.github.io/PyDCG/

DISK = 0
TRIANGLE = 0
assert(not DISK or not TRIANGLE)

gridsize = 2**int(argv[1]) # the first parameter specifies the gridsize (power of two)
k = int(argv[2]) # the second parameter specifies the "k", i.e., the size of the holes to count
n = int(argv[3]) # the third parameter specifies the number of points n 
n2 = n*n
num_tests = int(argv[4]) # the forth parameter specifies the number of samples

shape = argv[5] # fifth parameter specifies shape
if shape == "t": shape = "triangle"
if shape == "s": shape = "square"
if shape == "d": shape = "disk"
assert(shape in ["triangle","square","disk"])

N = range(n)

samples = []

for t in range(num_tests):
    while True:
        pts = []
        while len(pts) < n:
            x = randint(-gridsize,+gridsize)
            y = randint(-gridsize,+gridsize)
            if shape == "triangle":
                if x*x + y*y <= gridsize*gridsize: pts.append([x,y])
            elif shape == "disk": 
                if x <= y: pts.append([x,y])
            elif shape == "square": 
                pts.append([x,y])
            else:
                exit("shape not implemented")
           
        break

    # count empty k-holes
    if k > 3:
        empty_count = PyDCG.holesCpp.count_convex_rholes(pts,k)
    else:
        empty_count = PyDCG.holesCpp.countEmptyTriangs(pts)

    val = empty_count/float(n2)
    print t,"->",round(val,4)
    samples.append(val)

print 30*"-"+"summmary"+30*"-"

print "n:",n
print "k:",k
print "gridsize:",gridsize
print "shape:",shape

m = sum(samples)/num_tests
print "mean:\t",round(m,4)
s2 = sum((x-m)*(x-m) for x in samples)/(num_tests-1)
s = sqrt(s2)
print "s   :\t",round(s,4)


 
