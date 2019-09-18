"""
Copyright (c) 2006, Russell A. Paielli <http://RussP.us>
All rights reserved.

This file substitutes for scalar.py when the scalar class is switched
off for efficiency. It replaces the _scalar class with a built-in
numeric type such as float or int.

A user guide is available at <http://RussP.us/scalar.htm>. Send
comments, corrections, or suggestions to Russ.Paielli@gmail.com or use
the contact link at http://RussP.us.
"""

from __future__ import division # let 5/2 = 2.5 rather than 2 (caused problems)

from mathx import *

_convertTo={} # unit conversion factors

def _num(x):
    return x


def unit(name="", derived=1):

    _convertTo[name] = 1 / derived

    return derived

def unit_(name, equiv):
    "create a derived physical unit"

    _convertTo[name] = 1 / equiv

    return equiv

def AsUnit(arg, unit, fmt=""):

#    isarray = False
#    
#    try:
#        if len(arg) > 1: isarray = True
#    except:
#        pass
#    
#    if isarray:
#        string = "["
#        for x in arg: string += AsUnit(x, unit, fmt) + ", "
#        if arg: string = string[:-2]
#        string += "]"
#        return string

    temp = arg * _convertTo[unit]

    if fmt: out = fmt % temp
    else: out = str(temp)

    return out + " " + unit if unit else out

def unit_type(unit, utype=""):

    if not utype: return "unit unknown"
