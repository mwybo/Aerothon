#!/usr/bin/env python

"This file extends the standard math module with a few useful functions."

from __future__ import division # let 5/2 = 2.5 rather than 2 (caused problems)

from math import *

twoPi = 2 * pi
halfPi = pi / 2

def istr(x):
    "modified str (string) function, prints integer floats as ints"

    flt = float(str(x))
    if flt == int(flt): return str(int(flt))
    return str(x)

def normAngle(x):
    "add/subtract 2*pi to get angle in [-pi,+pi]"

    if x < -pi: return x + twoPi
    if x >  pi: return x - twoPi

    return x

def midAngle(a1, a2):
    "return angle halfway between a1 and a2"

    return a1 + normAngle(a2 - a1) / 2.

def sign(x):

    if x < 0: return -1
    if x > 0: return  1

    return 0

def nearmult(x, y=1): return int(x / y + (.5 if x > 0 else -.5)) * y # nearest multiple
def nearmultdif(x, y): return x - nearmult(x, y) # difference from nearest multiple
def nearmultdev(x, y): return abs(x - nearmult(x, y)) # deviation from nearest multiple
