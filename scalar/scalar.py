#!/usr/bin/env python

"""
Copyright (c) 2006, Russell A. Paielli <http://RussP.us>
All rights reserved.

See license file for terms of license.

This script implements scalar physical quantities with units. A user
guide is available at <http://RussP.us/scalar.htm>. Send comments,
corrections, or suggestions to Russ.Paielli@gmail.com or use the contact
link at http://RussP.us.
"""

from __future__ import division # let 5/2 = 2.5 rather than 2

from mathx import *
import numpy as npy

_convertTo = {} # dictionary for output conversions

class _scalar(object):
    "scalar physical quantities with units"

    #
    # Make _scalar dominant when multiplying with a numpy.array for left or right
    #
    __array_priority__ = 20.0
   
    class InconsistentUnits(Exception): pass
    class BadUnitName(Exception): pass # non-alphabetic character in unit name

    def __init__(self, num=1, units={}): # constructor

        self.num = num # number: float or int coefficient

        self.units = {} # storage for units

        for u in units: self.units[u] = units[u]

        self.cleanup()

    def __nonzero__(self): # for "truth-value" testing (conversion to boolean)

        if isinstance(self.num, npy.ndarray): return 1 if npy.all(self.num == 0) else 0
        if self.num: return 1
        return 0
    
    def __float__(self): # conversion to float (makes trig functions work)

        if not self.units: return float(self.num)
        if self.num == 0: return 0.

        raise _scalar.InconsistentUnits, \
            "can't convert " + str(self) + " to float"

    def __int__(self): # conversion to int

        if not self.units: return int(self.num)
        if self.num == 0: return 0

        raise _scalar.InconsistentUnits, \
            "can't convert " + str(self) + " to int"

    def __str__(self, star=0): # conversion to string for "print" output

        for unit in self.units:
            if not unit.isalpha(): raise _scalar.BadUnitName, \
               "non-alphabetic character in unit name: " + unit

        string = str(self.num) + self.printunits(star)
#        if not _is_zero(self.num): string += self.printunits(star)

        return string

    def __repr__(self): # string "representation" (can be executed to reconstruct scalar)

        return self.__str__(" * ")

    def copy(self): return _scalar(self.num, self.units)

    def checkunits(self, other):
        "assert consistent units for operations"

        self.cleanup()
        if isinstance(other, _scalar): other.cleanup()

        if _is_array_zero(self.num) or _is_array_zero(other): return # zero can have any units!
        if isinstance(other, npy.ndarray):
            sameunits = True
            for o in other:
                sameunits = sameunits and self.units == _units(o)
            if sameunits: return
        else:
            if self.units == _units(other): return

        raise _scalar.InconsistentUnits, str(self) + ", " + str(other)

    def cleanup(self):
        "remove units with zero exponent or no name"

        units = self.units

        for unit in units.keys(): # .keys() is needed to delete!
            if not unit or not units[unit]:
                del units[unit]

    def printunits(self, star=0):
        "construct unit string for output"

        if star: out = " * " # for repr output
        else: out = " " # for print output

        sep = "" # separator character between units
        units = self.units
        if not units: return ""

        ndenom = 0 # count of units in denominator (negative exponents)

        for unit in units:
            exp = units[unit]
            if exp < 0: ndenom += 1; continue # numerator first
            out += sep + unit
            if exp != 1: out += "**" + str(exp)
            sep = "*"

        if ndenom == 0: return out # no denominator
        if ndenom == 1: sep = "/" # one unit in denominator
        else: sep = "/(" # multiple units in denominator

        for unit in units:
            exp = units[unit]
            if exp >= 0: continue # denominator follows numerator
            out += sep + unit
            if exp != -1: out += "**" + str(-exp)
            sep = "*"

        if ndenom > 1: out += ")"

        return out

    # -- Methods to put any indexable type as scalar's num -----------

    def __getitem__(self,index):
        ''' returns a _scalar having the same unit as self,
            with a value being self's value sliced to index
        '''
        if hasattr(self.num, '__len__'):
            return _scalar(self.num[index], self.units)
        elif index == 0 or index == -1:
            return _scalar(self.num,self.units)
        else:
            raise IndexError

    def __setitem__(self,index,other):
        ''' makes a slice assignment on self's value based on index
            from value converted to self's unit
            raises DimensionError exception if self and value have incompatible units
        '''
        self.checkunits(other)
        if hasattr(self.num, '__len__'):
           self.num[index] = _num(other)
        elif index == 0 or index == -1:
           self.num = _num(other)
        else:
           raise IndexError

    def __len__(self):
        ''' returns the length of self's value
        '''
        if hasattr(self.num, '__len__'):
            return len(self.num)
        else:
            return NotImplemented
       
    def __neg__(self): # unary -
        return _scalar(-self.num, self.units)

    def __abs__(self): # absolute value
        return _scalar(abs(self.num), self.units)

    def __add__(self, other): # +

        if _is_zero(self): return _copy(other)

        self.checkunits(other)
        return _scalar(self.num + _num(other), self.units)

    def __sub__(self, other): # -

        if _is_zero(self): return -_copy(other)

        self.checkunits(other)
        return _scalar(self.num - _num(other), self.units)

    def __mul__(self, other): # *

#        if hasattr(other, "__len__"): return other.__mul__(self) # for lists

        result = _scalar(self.num * _num(other), self.units) # multiplication

        if isinstance(other,_scalar):
            units = result.units
            for unit, exp in other.units.items():
                if unit in units: units[unit] += exp
                else: units[unit] = exp

        result.cleanup()
        return _clean(result)

    def __div__(self, other):     # /
        return self.__truediv__(other)
    
    def __truediv__(self, other): # /

        result = _scalar(self.num / _num(other), self.units)
 
        if isinstance(other,_scalar):
            units = result.units
            for unit, exp in other.units.items():
                if unit in units: units[unit] -= exp
                else: units[unit] = -exp

        result.cleanup()
        return _clean(result)

    def __mod__(self, other): # %

        result.checkunits(other)
        return _scalar(self.num % _num(other), self.units)

    def __pow__(self, exp): # **

        result = _scalar(self.num**exp, self.units)
        if exp == 1: return result
        if exp == int(exp): exp = int(exp)
        for unit in result.units: result.units[unit] *= exp
        return result

    def __iadd__(self, other): return self + other;               # +=
    def __isub__(self, other): return self - other;               # -=
    def __imul__(self, other): return _clean(self * other);   # *=
    def __imod__(self, other): return self % other;            # %=
    def __ipow__(self, exp  ): return self**exp;             # **=

    def __idiv__ (self, other): return _clean(self / other);  # /=
    def __itruediv__ (self, other): return _clean(self / other);  # /=

    def __radd__(self, other): return self + other
    def __rmul__(self, other): return self * other

    def __rsub__(self, other): return -self + other
    def __rdiv__(self, other): return self**(-1) * other
    def __rtruediv__(self, other): return self**(-1) * other

    def __cmp__(self, other): # compare magnitude

        self.checkunits(other)

        if self.num > _num(other): return  1
        if self.num < _num(other): return -1

        return 0

    def _utuple(self):
        "construct hashable tuple from units"

        utup = self.units.items()
        utup.sort()
        return tuple(utup)

"""The following functions would normally be member functions, but public
member functions are not allowed for the scalar class because all
functions must work on built-in types also (when the scalar class is
turned off for efficiency)."""

def _clean(x):

    if isinstance(x, _scalar): 
        return x.num if not x.units else x        
#        return x.num if not x.units or _is_zero(x) else x        
    return x

def _num(x):

    if isinstance(x, _scalar): return x.num
    if isinstance(x, npy.ndarray): 
        try:
            return npy.array([xi.num for xi in x])
        except:
            return x
        
    return x

def _is_zero(x):
    ''' Determines if x is zero, returns false if x or x.num is an array'''
    if isinstance(x, _scalar): 
        if isinstance(x.num, npy.ndarray):
            return False
        else:
            return x.num == 0
        
    if isinstance(x, npy.ndarray): return False
    return x == 0

def _is_array_zero(x):
    ''' Determines if x is zero, and returns true of an array is all zeros'''
    if isinstance(x, _scalar): 
        if isinstance(x.num, npy.ndarray):
            return npy.all(x.num == 0)
        else:
            return x.num == 0
        
    if isinstance(x, npy.ndarray): return npy.all(x == 0)
    return x == 0

def _units(x):

    if not isinstance(x, _scalar): return {}
    units = {}
    for u in x.units: units[u] = x.units[u]
    return units

def _copy(x):

    if isinstance(x, _scalar): return x.copy()
    return x

def AsUnit(arg, unit="", fmt=""):
    "specify output unit and numeric AsUnit"

    #if hasattr (arg, "__len__"): # for list types

    isarray = False
    
    try:
        if len(arg) > 1: isarray = True
    except:
        pass
    
    if isarray:
        string = "["
        for x in arg: string += AsUnit(x, unit, fmt) + ", "
        if arg: string = string[:-2]
        string += "]"
        return string

    temp = arg * _convertTo[unit] if unit else arg
    num = _num(temp)

    if fmt: out = fmt % num
    else: out = str(num)

    if isinstance (temp, _scalar):

        for unit0 in _units(arg):
            if not unit0.isalpha(): raise _scalar.BadUnitName, \
               "non-alphabetic character in unit name: " + unit0

        if temp.printunits().strip() != unit and temp:
            raise _scalar.InconsistentUnits, str(temp) + ", " + unit

    return out + " " + unit if unit else out

def sqrt(x): return x**0.5 # sqrt for scalar and float

def hypot(x, y): return (x**2 + y**2)**0.5 # hypot for scalar and float

_atan2 = atan2 # save reference to original atan2 function

def atan2(y, x): # atan2 for scalar and float

    ax = abs(x)
    ay = abs(y)

    if ax > ay: return _atan2(y/ax, sign(x))
    return _atan2(sign(y), x/ay)

base_units = [] # storage for reference list of base units
derived_units = [] # storage for reference list of derived units

def unit(name="", equiv=0):
    "create a base physical unit"

    if equiv: return unit_(name, equiv)

    base_units.append(name)

    _convertTo[name] = 1

    return _scalar(1, {name: 1})

def unit_(name, equiv):
    "create a derived physical unit"

    derived_units.append(name)

    _convertTo[name] = _scalar(1, {name: 1}) / equiv

    return equiv

def unit_type(unit, utype="", qdict={}):
    "set or get the physical type of a unit (time, length, etc.)"

    utuple = unit._utuple() # hashable tuple representation of units

    if utype: qdict[utuple] = utype; return

    if utuple not in qdict: return "unknown"
    return qdict[utuple]
