#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 10 17:12:27 2018

@author: egissi
"""

import math

EPSILON = 1e-07
EPSILON2 = 1e-05


class Vector(object):
    def __init__(self, *args):  # FIXME simplify
        self.x, self.y, self.z = 0., 0., 0.
        if len(args) == 3:  # Vector(1,2,3)
            self.x, self.y, self.z = args[0], args[1], args[2]
        elif len(args) == 1:  # Vector([1,2,3])
            a = args[0]
            if isinstance(a, dict):
                self.x, self.y, self.z = \
                    a.get('x', 0.0), a.get('y', 0.0), a.get('z', 0.0)
            elif a is not None and len(a) == 3:
                self.x, self.y, self.z = a[0], a[1], a[2]

    def __repr__(self):
        return 'Vector({0:.3f}, {1:.3f}, {2:.3f})'.format(
                self.x, self.y, self.z
                )

    def clone(self):
        return Vector(self.x, self.y, self.z)

    def negated(self):
        return Vector(-self.x, -self.y, -self.z)

    def __neg__(self):
        return self.negated()

    def plus(self, a):
        return Vector(self.x+a.x, self.y+a.y, self.z+a.z)

    def __add__(self, a):
        return self.plus(a)

    def minus(self, a):
        return Vector(self.x-a.x, self.y-a.y, self.z-a.z)

    def __sub__(self, a):
        return self.minus(a)

    def times(self, a):
        return Vector(self.x*a, self.y*a, self.z*a)

    def __mul__(self, a):
        return self.times(a)

    def dividedBy(self, a):
        return Vector(self.x/a, self.y/a, self.z/a)

    def __truediv__(self, a):
        return self.dividedBy(float(a))

    def __div__(self, a):
        return self.dividedBy(float(a))

    def dot(self, a):
        return self.x*a.x + self.y*a.y + self.z*a.z

    def lerp(self, a, t):
        "Linear interpolation from self to a"
        return self.plus(a.minus(self).times(t))

    def length(self):
        return math.sqrt(self.dot(self))

    def unit(self):
        return self.dividedBy(self.length())

    def cross(self, a):
        return Vector(
            self.y * a.z - self.z * a.y,
            self.z * a.x - self.x * a.z,
            self.x * a.y - self.y * a.x,
            )

    def __getitem__(self, key):
        return (self.x, self.y, self.z)[key]

    def __setitem__(self, key, value):
        ll = [self.x, self.y, self.z]
        ll[key] = value
        self.x, self.y, self.z = ll

    def __len__(self):
        return 3

    def __iter__(self):
        return iter((self.x, self.y, self.z))

    def __eq__(self, other):
        return \
            math.isclose(self.x, other.x, rel_tol=EPSILON) and \
            math.isclose(self.y, other.y, rel_tol=EPSILON) and \
            math.isclose(self.z, other.z, rel_tol=EPSILON)

    def __hash__(self):
        return hash((self.x, self.y, self.z))

    def isZero(self):
        return abs(self.x) < EPSILON and \
               abs(self.y) < EPSILON and \
               abs(self.z) < EPSILON

    def is_collinear(self, b, c):
        "Return true if a, b, and c all lie on the same line."
        # Proposed elsewhere, to correctly use epsilon FIXME FIXME
        # area_tri = abs((c-self).cross(b-self).length())
        # area_square = max((c-self).length() ** 2,
        #   (b-self).length() ** 2, (c-b).length() ** 2)
        # if 2* area_tri < EPSILON * area_square: return True
        return abs((c-self).cross(b-self).length()) < EPSILON

    def is_within(self, p, r):
        "Return true if q is between p and r (inclusive)."
        return (p.x <= self.x <= r.x or r.x <= self.x <= p.x) and \
               (p.y <= self.y <= r.y or r.y <= self.y <= p.y) and \
               (p.z <= self.z <= r.z or r.z <= self.z <= p.z)

if __name__ == "__main__":
    import doctest
    doctest.testmod()
