#!/usr/bin/env python

from __future__ import division # let 5/2 = 2.5 rather than 2

from math import *

_halfPi = pi / 2

class vector(object):
    "basic vector operations"

    class InconsistentLengths(Exception): pass

    def __init__(self, vec): self.vec = list(vec)

    def __len__(self): return len(self.vec)
    def __getitem__(self, i): return self.vec[i]
    def __setitem__(self, i, x): self.vec[i] = x
    def __delitem__(self, i): del self.vec[i]
    def __getslice__(self, i, j): return vector(self.vec[i:j])
    def __setslice__(self, i, j, s): self.vec[i:j] = s
    def __delslice__(self, i, j): del self.vec[i:j]
    def __contains__(self, obj): return obj in self.vec

    def __checklength (self, other):

        if len(self) == len(other): return
        strx = str(len(self)) + ", " + str(len(other))
        raise vector.InconsistentLengths, strx

    def __nonzero__(self): # for "truth-value" testing (conversion to boolean)

        if self.vec: return 1
        return 0

    def __str__(self): # convert to string for "print" statement

        out = "["
        for x in self: out += str(x) + ", "
        if self.vec: out = out[:-2]
        out += "]"
        return out

    def __repr__(self): return "vector(" + str(self) + ")"

    def mag(self): # vector magnitude (length)

        Sum = 0.
        for x in self: Sum += x**2
        return Sum**0.5

    def append(self, item): self.vec.append(item)

    def atan2(self): return atan2(self[1], self[0])

    def AsUnit(self, fmt): # convert output

        string = "["
        for x in range(len(self)): string += fmt % x + ", "
        if self.vec: string = string[:-2]
        string += "]"
        return string

    def __neg__(self): # unary -

        result = self[:]
        for i in range(len(result)): result[i] = -result[i]
        return result

    def __iadd__(self, other): # +=

        self.__checklength(other)
        for i in range(len(self)): self[i] += other[i]
        return self

    def __isub__(self, other): # -=

        self.__checklength(other)
        for i in range(len(self)): self[i] -= other[i]
        return self

    def __imul__(self, other): # *=

        for i in range(len(self)): self[i] *= other
        return self

    def __itruediv__(self, other): # /=

        for i in range(len(self)): self[i] /= other
        return self

    def __imod__(self, other): # %=

        for i in range(len(self)): self[i] %= other
        return self

    def __add__(self, other): # +

        result = self[:]
        result += other
        return result

    def __sub__(self, other): # -

        result = self[:]
        result -= other
        return result

    def __mul__(self, other): # *

        result = self[:]
        result *= other
        return result

    def __mod__(self, other): # %

        result = self[:]
        result %= other
        return result

    def __truediv__ (self, other): # /

        result = self[:]
        result /= other
        return result

    def __radd__(self, other): return self + other
    def __rmul__(self, other): return self * other
    def __rsub__(self, other): return other - self
    def __rtruediv__(self, other): return other / self

    def dotProduct (self, vec2): # vector dot product

        self.__checklength(vec2)
        prod = 0.
        for i in range(len(self)): prod += self[i] * vec2[i]
        return prod

def isZero(vec):

    for i in range(len(vec)):
        if vec[i]: return 0
    return 1

def isNonZero(vec):

    for i in range(len(vec)):
        if vec[i]: return 1
    return 0

def sumSquare(vec):

    result = 0
    for i in range(len(vec)): result += vec[i]**2
    return result

def magnitude(vec):

    result = 0
    for i in range(len(vec)): result += vec[i]**2
    return result**0.5

def dotProduct(vec1, vec2):

    result = 0
    for i in range(len(vec1)): result += vec1[i] * vec2[i]
    return result

def negative(vec):

    result = vec[:]
    for i in range(len(vec)): result[i] = -result[i]
    return result

def add(vec1, vec2):

    result = vec1[:]
    for i in range(len(vec1)): result[i] += vec2[i]
    return result

def subtract(vec1, vec2):

    result = vec1[:]
    for i in range(len(vec1)): result[i] -= vec2[i]
    return result

def multiply(vec, scalar):

    result = vec[:]
    for i in range(len(vec)): result[i] *= scalar
    return result

def divide(vec, scalar):

    result = vec[:]
    for i in range(len(vec)): result[i] /= scalar
    return result
