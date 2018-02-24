#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 10 17:12:27 2018

@author: egissi
"""

import math
import array
import textwrap

EPSILON = 1e-4  # FIXME different EPSILON for different applications
EPSILON_CUT = 1e-04


class Vector(object):
    def __init__(self, *args):
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
        """
        Linear interpolation from self to a.
        >>> Vector(0,0,1).lerp(Vector(10,0,1),.7)
        Vector(7.000, 0.000, 1.000)
        """
        return self.plus(a.minus(self).times(t))

    def length(self):
        return math.sqrt(self.dot(self))

    def squarelength(self):
        return self.dot(self)

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

    def __hash__(self):
        return hash((self.x, self.y, self.z))

    def __eq__(self, other):
        return \
            math.isclose(self.x, other.x, rel_tol=EPSILON) and \
            math.isclose(self.y, other.y, rel_tol=EPSILON) and \
            math.isclose(self.z, other.z, rel_tol=EPSILON)

    def is_zero(self, multiplier=1.):  # FIXME Epsilon
        # abs(a-b) <= max( rel_tol * max(abs(a), abs(b)), abs_tol )
        return \
            math.isclose(self.x, 0., rel_tol=EPSILON * multiplier) and \
            math.isclose(self.y, 0., rel_tol=EPSILON * multiplier) and \
            math.isclose(self.z, 0., rel_tol=EPSILON * multiplier)

    def is_collinear(self, b, c, multiplier=1.):
        """
        Return true if self, b, and c all lie on the same line.
        >>> a, b, c = Vector(0,0,0), Vector(1,0,0), Vector(2,0,0)
        >>> a.is_collinear(b, c)  # Collinear
        True
        >>> a, b, c = Vector(0,0,0), Vector(0,-1,0), Vector(0,-2,1)
        >>> a.is_collinear(b, c)  # Not collinear
        False
        >>> a, b, c = Vector(0,0,0), Vector(0,0,0), Vector(0,0,0)
        >>> a.is_collinear(b, c)  # Same point
        True
        >>> a, b, c = Vector(0,0,0), Vector(1,0,0), Vector(2,EPSILON,0)
        >>> a.is_collinear(b, c)
        True
        """
        # Simpler
        # return b.minus(self).cross(c.minus(self)).is_zero()
        # Proposed elsewhere, to correctly use epsilon FIXME
        c_s = c.minus(self)
        b_s = b.minus(self)
        if c_s.is_zero() or b_s.is_zero():
            return True
        area_tri2 = c_s.cross(b_s).length()  # FIXME go to squarelenght()
        area_square = max(
                c_s.length() ** 2,
                b_s.length() ** 2,
                (c.minus(b)).length() ** 2
                )
        if area_tri2 < EPSILON * area_square * multiplier:
            return True
        return False

    def is_within(self, p, r):
        """
        Return true if q is between p and r (inclusive).
        >>> Vector(0,0,0).is_within(Vector(1,0,0), Vector(2,0,0))
        False
        >>> Vector(1.5,0,0).is_within(Vector(1,0,0), Vector(2,0,0))
        True
        """
        return (p.x <= self.x <= r.x or r.x <= self.x <= p.x) and \
               (p.y <= self.y <= r.y or r.y <= self.y <= p.y) and \
               (p.z <= self.z <= r.z or r.z <= self.z <= p.z)

    def is_within_tri(self, a, b, c, normal=None):
        """
        Test if a point is in a triangle.
        A pre-calculated normal can be sent
        >>> Vector(1,1,0).is_within_tri(\
                Vector(0,0,0),Vector(2,0,0),Vector(0,2,0))
        True
        >>> Vector(2,0,0).is_within_tri(\
                Vector(0,0,0),Vector(2,0,0),Vector(0,2,0))
        True
        >>> Vector(3,0,0).is_within_tri(\
                Vector(0,0,0),Vector(2,0,0),Vector(0,2,0))
        False
        >>> Vector(3,1,0).is_within_tri(\
                Vector(0,3,0),Vector(0,0,0),Vector(3,0,0))
        False
        """
        # Test tri bounding box
        p = Vector(
                max(a.x, b.x, c.x),
                max(a.y, b.y, c.y),
                max(a.z, b.z, c.z),
                )
        q = Vector(
                min(a.x, b.x, c.x),
                min(a.y, b.y, c.y),
                min(a.z, b.z, c.z),
                )
        if not self.is_within(p, q):
            return False
        # Test cross products
        if normal is None:
            normal = Plane.from_points((a, b, c)).normal
        if all((
                b.minus(a).cross(self.minus(a)).dot(normal) >= 0.,
                c.minus(b).cross(self.minus(b)).dot(normal) >= 0.,
                a.minus(c).cross(self.minus(c)).dot(normal) >= 0.,
                )):  # FIXME EPSILON precision!
            return True
        else:
            return False


class Plane():
    def __init__(self, normal, distance):
        self.normal = Vector(normal).unit()
        self.distance = distance  # distance from (0, 0, 0)

    def __repr__(self):
        return 'Plane(normal={}, distance={:.3f})'.format(
                self.normal, self.distance
                )

    def clone(self):
        return Plane(self.normal.clone(), self.distance)

    @classmethod
    def from_points(cls, points):
        """
        Get a Plane from a list of points, after checking for collinearity.
        Robust algorithm, good for slightly concave polygons.
        If collinear, return a None.
        >>> Plane.from_points(((0,0,5,),(1,0,5),(2,0,5),(0,1,5)))   # Convex
        Plane(normal=Vector(0.000, 0.000, 1.000), distance=5.000)
        >>> Plane.from_points(((0,0,5,),(1,.1,5),(2,0,5),(0,1,5)))  # Concave
        Plane(normal=Vector(0.000, 0.000, 1.000), distance=5.000)
        >>> Plane.from_points(((1,0,0,),(0,1,0),(1,0,1),))          # Inclined
        Plane(normal=Vector(0.707, 0.707, 0.000), distance=0.707)
        >>> Plane.from_points(((0,0,0,),(1,0,0),(2,0,0),(3,0,0))) \
            # Collinear # doctest: +NORMALIZE_WHITESPACE
        Traceback (most recent call last):
        ...
        Exception: ('Could not find a plane, points:',
                    ((0, 0, 0), (1, 0, 0), (2, 0, 0), (3, 0, 0)))
        """
        len_points = len(points)
        tot_normal = Vector()
        for i, a in enumerate(points):
            a = Vector(a)
            b = Vector(points[(i + 1) % len_points])
            c = Vector(points[(i + 2) % len_points])
            tot_normal += b.minus(a).cross(c.minus(a))
        if tot_normal.is_zero():
            raise Exception('Could not find a plane, points:', points)
        normal = tot_normal.unit()
        return Plane(normal, a.dot(normal))

    def flip(self):
        """
        Flip self normal.
        >>> p = Plane(normal=(1,0,0), distance=5); p.flip() ; p
        Plane(normal=Vector(-1.000, -0.000, -0.000), distance=-5.000)
        """
        self.normal = -self.normal
        self.distance = -self.distance


class Geom():
    """
    Representation of a polygonal geometry
    """
    def __init__(self, verts=None, polygons=None, surfids=None, hid=None):
        self.hid = hid  # Geom name
        if verts is None:
            verts = ()
        if polygons is None:
            polygons = ()
        try:
            # Array of vertices coordinates, eg. [1.,2.,3., 2.,3.,4., ...]
            self.verts = array.array('f', verts)
            # List of list of polygon connectivities, eg. [[0,1,2,4], ...]
            self.polygons = [list(p[:]) for p in polygons]
            # List of polygon surfid indexes, eg. [0,0,1,2,5, ...]
        except TypeError:
            raise Exception('Bad Geom(), hid:', hid)
        # Set surfids
        if surfids is None:
            self.surfids = ['White' for polygon in self.polygons]
        else:
            self.surfids = list(surfids)
        if len(self.surfids) != len(self.polygons):
            raise Exception('Bad surfids in Geom(), hid:', hid)
        # Set normals
        self.normals = []
        if self.polygons and self.verts:
            self.update_normals()

    def __repr__(self):
        strverts = ('{:.3f}'.format(v) for v in self.verts)
        strverts = list(zip(*[iter(strverts)] * 3))  # Join: ((1.,2.,3.,), ...)
        strverts = ',  '.join((','.join(v) for v in strverts))
        return 'Geom(\n    ({}),\n    {},\n    )'.format(
                strverts, self.polygons,
                )

    def clone(self):
        geom = Geom(
                verts=self.verts[:],
                polygons=[p[:] for p in self.polygons],
                surfids=self.surfids[:],
                hid=self.hid,
                )
        # Never recalc normals, if possible
        geom.normals = self.normals[:]
        return geom

    # Modify geom

    def append(self, geom):
        """
        Append other geom to self.
        >>> g = Geom((-1,-1,0, 1,-1,0, 0,1,0, 0,0,1), \
                     ((2,1,0), (0,1,3), (1,2,3), (2,0,3)) )  # Good tet
        >>> h = Geom((-1,-1,2, 1,-1,2, 0,1,2, 0,0,1), \
                     ((0,1,2), (3,1,0), (3,2,1), (3,0,2)) )  # but upside-down
        >>> g.append(h); g  # doctest: +NORMALIZE_WHITESPACE
        Dup verts removed: 1
        [0, 1, 2, 3, 4, 5, 6, 7]
        Geom(
            (-1.000,-1.000,0.000,  1.000,-1.000,0.000,  0.000,1.000,0.000,
            0.000,0.000,1.000,  -1.000,-1.000,2.000,  1.000,-1.000,2.000,
            0.000,1.000,2.000),
            [[2, 1, 0], [0, 1, 3], [1, 2, 3], [2, 0, 3], [4, 5, 6],
            [3, 5, 4], [3, 6, 5], [3, 4, 6]],
            )
        """
        # Join cut border of self and geom
        border0_halfedges = self.get_border_halfedges()
        iverts0 = [bh[0] for bh in border0_halfedges]
        border1_halfedges = geom.get_border_halfedges()
        iverts1 = [bh[0] for bh in border1_halfedges]
        # Merge iverts0 and ivert1
        while iverts0:
            ivert0 = iverts0.pop()
            vert0 = self.get_vert(ivert0)
            for ivert1 in iverts1:
                vert1 = geom.get_vert(ivert1)
                if vert0.minus(vert1).length() < EPSILON_CUT * 2:  # FIXME go to squarelenght
                    geom.update_vert(ivert1, vert0)
                    iverts1.remove(ivert1)
                    break
        # Check if all went well
        if iverts0:
            print("Unmerged iverts0:", iverts0)
        if iverts1:
            print("Unmerged iverts1:", iverts1)
        # Extend verts of self with geom's
        original_nverts = self.get_nverts()
        self.verts.extend(geom.verts)
        # Extend polygons of self with geom's
        original_npolygons = self.get_npolygons()
        self.polygons.extend(geom.polygons)
        # Relink geom polygons to new iverts
        for i, polygon in enumerate(self.polygons[original_npolygons:]):
            for j, _ in enumerate(polygon):
                polygon[j] += original_nverts
        # Extend self surfids with geom's
        self.surfids.extend(geom.surfids)
        # Extend self normals with geom's
        self.normals.extend(geom.normals)
        # Merge duplicate verts
        self.merge_duplicated_verts()
        return self.get_ipolygons()  # FIXME why?

    def flip(self):
        """
        Flip all polygon normals.
        >>> g = Geom((0,0,1, 1,0,1, 2,0,1, 0,1,1,), ((0,1,2,3), ))
        >>> g.flip(); g.get_polygon_normal(0); g.get_polygon(0)
        Vector(-0.000, -0.000, -1.000)
        [3, 2, 1, 0]
        """
        for polygon in self.polygons:
            polygon.reverse()
        for ipolygon, normal in enumerate(self.normals):
            self.normals[ipolygon] = normal.negated()

    def update_normals(self):
        """
        Update polygon normals
        """
        self.normals = [Plane.from_points(
                self.get_polygon_verts(ipolygon)
                ).normal for ipolygon, polygon in enumerate(self.polygons)]

    # Numbers

    def get_npolygons(self):
        """
        Get the len of polygons
        >>> g = Geom((-1,-1,0, 1,-1,0, 0,1,0, 0,0,1), \
                     ((2,1,0), (0,1,3), (1,2,3), (2,0,3)) )  # Good tet
        >>> g.get_npolygons()
        4
        """
        return len(self.polygons)

    def get_ipolygons(self):
        """
        Get the list of ipolygon
        >>> g = Geom((-1,-1,0, 1,-1,0, 0,1,0, 0,0,1), \
                     ((2,1,0), (0,1,3), (1,2,3), (2,0,3)) )  # Good tet
        >>> g.get_ipolygons()
        [0, 1, 2, 3]
        """
        return list(range(len(self.polygons)))

    # Get polygon

    def get_polygon(self, ipolygon):
        """
        Get ipolygon connectivity.
        >>> g = Geom((0,0,1, 1,0,1, 2,0,1, 0,1,1,), ((0,1,2,3), ))
        >>> g.get_polygon(0)
        [0, 1, 2, 3]
        """
        return self.polygons[ipolygon]

    def get_polygon_verts(self, ipolygon):
        """
        Get ipolygon verts.
        >>> g = Geom((-1,-1,0, 1,-1,0, 0,1,0, 0,0,1), \
                     ((2,1,0), (0,1,3), (1,2,3), (2,0,3)) )  # Good tet
        >>> g.get_polygon_verts(0)  # doctest: +NORMALIZE_WHITESPACE
        [Vector(0.000, 1.000, 0.000), Vector(1.000, -1.000, 0.000),
         Vector(-1.000, -1.000, 0.000)]
        """
        return [self.get_vert(ivert) for ivert in self.get_polygon(ipolygon)]

    def get_polygon_surfid(self, ipolygon):
        """
        Get ipolygon surfid.
        >>> g = Geom((-1,-1,0, 1,-1,0, 0,1,0, 0,0,1), \
                     ((2,1,0), (0,1,3), (1,2,3), (2,0,3)), \
                     (0,3,2,1), )  # Good tet
        >>> g.get_polygon_surfid(1)
        3
        """
        return self.surfids[ipolygon]

    def get_polygon_normal(self, ipolygon):
        """
        Get normal of ipolygon.
        >>> g = Geom((0,0,1, 1,0,1, 2,0,1, 0,1,1,), ((0,1,2,3), ))
        >>> g.get_polygon_normal(0)
        Vector(0.000, 0.000, 1.000)
        """
        return self.normals[ipolygon]

    def get_polygon_plane(self, ipolygon):
        """
        Get plane containing ipolygon in a robust way.
        >>> g = Geom((0,0,1, 1,0,1, 2,0,1, 0,1,1,), ((0,1,2,3), ))
        >>> g.get_polygon_plane(0)
        Plane(normal=Vector(0.000, 0.000, 1.000), distance=1.000)
        """
        verts = self.get_polygon_verts(ipolygon)
        normal = self.normals[ipolygon]
        tot_distance = 0.
        for vert in verts:
            tot_distance += vert.dot(normal)
        distance = tot_distance / len(verts)
        return Plane(normal, distance)

    # Modify polygon

    def remove_polygon(self, ipolygon):
        """
        Remove a polygon. Change ipolygon references!
        """
        print("Polygon removed:", ipolygon)
        del self.polygons[ipolygon]
        del self.surfids[ipolygon]
        del self.normals[ipolygon]

    def update_polygon(self, ipolygon, polygon, normal=None):
        """
        Update ipolygon connectivity.
        >>> g = Geom((-1,-1,0, 1,-1,0, 0,1,0, 0,0,1), \
                     ((2,1,0), (0,1,3), (1,2,3), (2,0,3)) )  # Good tet
        >>> g.update_polygon(0, (0,1,2), (1,0,0)); g.get_polygon(0)
        0
        [0, 1, 2]
        >>> g.get_polygon_normal(0)
        Vector(1.000, 0.000, 0.000)
        """
        self.polygons[ipolygon] = list(polygon)
        if normal is not None:
            self.normals[ipolygon] = Vector(normal)
        return ipolygon

    def append_polygon(self, polygon, surfid, normal):
        """
        Append a polygon.
        >>> g = Geom((-1,-1,0, 1,-1,0, 0,1,0, 0,0,1), \
                     ((2,1,0), (0,1,3), (1,2,3), (2,0,3)) )  # Good tet
        >>> g.append_polygon((3,1,0), 1, (1,0,0)); g.get_polygon(4)
        4
        [3, 1, 0]
        >>> g.get_polygon_surfid(4); g.get_polygon_normal(4)
        1
        Vector(1.000, 0.000, 0.000)
        """
        self.polygons.append(list(polygon))
        self.surfids.append(surfid)
        self.normals.append(Vector(normal))
        ipolygon = self.get_npolygons() - 1
        return ipolygon

    def merge_polygons_to_concave(self, ipolygons):
        """
        Merge coplanar polygons (online) with same surfid to concave polygons.
        No check is performed on coplanarity. This is used by the bsp tree.
        """
        # Build surfid_to_polygons dict FIXME put in a def used twice (to_OBJ)
        surfid_to_ipolygons = {}
        for ipolygon in ipolygons:
            surfid = self.get_polygon_surfid(ipolygon)
            try:
                surfid_to_ipolygons[surfid].append(ipolygon)
            except KeyError:
                surfid_to_ipolygons[surfid] = [ipolygon, ]
        # For each surfid merge polygons
        for surfid, surfid_ipolygons in surfid_to_ipolygons.items():
            border_loops = self.get_border_loops(surfid_ipolygons)
            for loop, loop_ipolygons in border_loops.items():
                self.update_polygon(loop_ipolygons[0], loop)  # FIXME if not coplanar, set new normal
                for ipolygon in loop_ipolygons[1:]:
                    ipolygons.remove(ipolygon)

    def split_polygon(self, ipolygon, plane, coplanar_front,
                      coplanar_back, front, back):
        """
        Split ipolygon by a plane. Put the fragments in the inline lists.
        Add cut_ivert to bordering polygons.
        >>> g = Geom((-1,-1,0, 1,-1,0, 1,1,0, -1,1,0, -3, 1,0, -3,-1,0, \
                       3,-1,0, 3, 1,0, 1,3,0, -1,3,0, -1,-3,0,  1,-3,0),\
                     ((0,1,2,3), (5,0,3,4), (1,6,7,2), (3,2,8,9),\
                     (10,11,1,0)), (0,1,2,3,4), \
                    )  # Open clover on z=0, n=+k
        >>> h = g.clone()
        >>> #  y ↑
        >>> #   9─8
        >>> #   │ │
        >>> # 4─3─2─7
        >>> # │ │∙│ │→ x
        >>> # 5─0─1─6
        >>> #   │ │
        >>> #  10-11
        >>> coplanar_front, coplanar_back, front, back = [], [], [], []
        >>> g.split_polygon(0, Plane((1,0,0),0), \
                            coplanar_front, coplanar_back, front, back)
        >>> coplanar_front, coplanar_back, front, back
        ([], [], [0], [5])
        >>> g.to_OBJ('../test/clover.obj')
        to_OBJ: ../test/clover.obj
        >>> g  # doctest: +NORMALIZE_WHITESPACE
        Geom(
            (-1.000,-1.000,0.000,  1.000,-1.000,0.000,  1.000,1.000,0.000,
             -1.000,1.000,0.000,  -3.000,1.000,0.000,  -3.000,-1.000,0.000,
             3.000,-1.000,0.000,  3.000,1.000,0.000,  1.000,3.000,0.000,
             -1.000,3.000,0.000,  -1.000,-3.000,0.000,  1.000,-3.000,0.000,
             0.000,-1.000,0.000,  0.000,1.000,0.000),
            [[12, 1, 2, 13], [5, 0, 3, 4], [1, 6, 7, 2], [3, 13, 2, 8, 9],
            [10, 11, 1, 12, 0], [0, 12, 13, 3]],
            )
        >>> coplanar_front, coplanar_back, front, back = [], [], [], []
        >>> h.split_polygon(0, Plane((0,1,0),0), \
                            coplanar_front, coplanar_back, front, back)
        >>> coplanar_front, coplanar_back, front, back
        ([], [], [0], [5])
        >>> h  # doctest: +NORMALIZE_WHITESPACE
        Geom(
            (-1.000,-1.000,0.000,  1.000,-1.000,0.000,  1.000,1.000,0.000,
            -1.000,1.000,0.000,  -3.000,1.000,0.000,  -3.000,-1.000,0.000,
            3.000,-1.000,0.000,  3.000,1.000,0.000,  1.000,3.000,0.000,
            -1.000,3.000,0.000,  -1.000,-3.000,0.000,  1.000,-3.000,0.000,
            1.000,0.000,0.000,  -1.000,0.000,0.000),
            [[12, 2, 3, 13], [5, 0, 13, 3, 4], [1, 6, 7, 2, 12], [3, 2, 8, 9],
            [10, 11, 1, 0], [0, 1, 12, 13]],
            )
        >>> coplanar_front, coplanar_back, front, back = [], [], [], []
        >>> g.split_polygon(0, Plane((0,0,1),0), \
                            coplanar_front, coplanar_back, front, back)
        >>> coplanar_front, coplanar_back, front, back
        ([0], [], [], [])
        >>> coplanar_front, coplanar_back, front, back = [], [], [], []
        >>> g.split_polygon(0, Plane((0,0,-1),0), \
                            coplanar_front, coplanar_back, front, back)
        >>> coplanar_front, coplanar_back, front, back
        ([], [0], [], [])
        """
        # Init
        COPLANAR = 0  # vertex of polygon in EPSILON_CUT distance f plane
        FRONT = 1     # vertex of polygon in front of the plane
        BACK = 2      # vertex of polygon at the back of the plane
        SPANNING = 3  # spanning polygon

        polygon = self.get_polygon(ipolygon)
        polygon_type = 0
        ivert_types = []
        polygon_nverts = len(polygon)
        polygon_surfid = self.get_polygon_surfid(ipolygon)
        polygon_normal = self.get_polygon_normal(ipolygon)

        # The edges to be split,
        # eg. {(2,3): 1} with {(ivert0,ivert1): cut_ivert, ...}
        # The opposite of the split edge is sent for easier search
        spl_edges = {}

        # Calc the distance between the ivert and the splitting plane
        # then classify the ivert, and update classification of the face
        for ivert in polygon:
            # Classify ivert using vert-plane distance
            distance = plane.normal.dot(self.get_vert(ivert)) - plane.distance
            ivert_type = -1
            if distance < -EPSILON_CUT * abs(plane.distance):  # FIXME
                ivert_type = BACK
            elif distance > EPSILON_CUT * abs(plane.distance):  # FIXME
                ivert_type = FRONT
            else:
                ivert_type = COPLANAR
            # Register ivert classification
            ivert_types.append(ivert_type)
            # Update polygon classification
            polygon_type |= ivert_type

        # Put the polygon in the correct list
        if polygon_type == COPLANAR:
            # Same or opposite normal?
            if plane.normal.dot(polygon_normal) > 0:
                coplanar_front.append(ipolygon)
            else:
                coplanar_back.append(ipolygon)
        elif polygon_type == FRONT:
            front.append(ipolygon)
        elif polygon_type == BACK:
            back.append(ipolygon)
        elif polygon_type == SPANNING:
            front_iverts = []
            back_iverts = []
            for i, ivert0 in enumerate(polygon):
                # Get the edge ivert0-ivert1
                j = (i+1) % polygon_nverts
                ivert1 = polygon[j]
                ivert0_type = ivert_types[i]
                ivert1_type = ivert_types[j]
                # Put ivert0 in the right lists
                # to build the new edge
                if ivert0_type != BACK:
                    front_iverts.append(ivert0)
                if ivert0_type != FRONT:
                    if ivert0_type != BACK:
                        back_iverts.append(ivert0)
                    else:
                        back_iverts.append(ivert0)
                if (ivert0_type | ivert1_type) == SPANNING:
                    # The edge is spanning, calc new vert
                    vert0 = self.get_vert(ivert0)
                    vert1 = self.get_vert(ivert1)
                    t = (plane.distance - plane.normal.dot(vert0)) \
                        / plane.normal.dot(vert1 - vert0)
                    cut_vert = vert0.lerp(vert1, t)
                    cut_ivert = self.append_vert(cut_vert)
                    # Register the split for domino to bordering polygons
                    spl_edges[(ivert1, ivert0)] = cut_ivert  # opposite!
                    # Append the new_vert to the right list
                    front_iverts.append(cut_ivert)
                    back_iverts.append(cut_ivert)

            # Update and append new polygons
            updated = False
            if len(front_iverts) >= 3:
                updated = True
                new_ipolygon = self.update_polygon(ipolygon, front_iverts)
                front.append(new_ipolygon)

            if len(back_iverts) >= 3:
                if updated:
                    new_ipolygon = self.append_polygon(
                            back_iverts, polygon_surfid, polygon_normal,
                            )
                else:
                    updated = True
                    new_ipolygon = self.update_polygon(ipolygon, back_iverts)
                back.append(new_ipolygon)

            # Add cut_vert to bordering polygons
            halfedges = self.get_halfedges()
            for spl_edge, cut_ivert in spl_edges.items():
                # Get the bordering polygon that is split by cut_ivert
                spl_ipolygon = halfedges.get(spl_edge, None)
                if spl_ipolygon is None:  # there is a border
                    continue
                # Calc and update the bordering polygon
                spl_polygon = self.get_polygon(spl_ipolygon)
                i = spl_polygon.index(spl_edge[0])  # find right edge
                spl_polygon.insert(i+1, cut_ivert)  # inject cut_ivert
                self.update_polygon(spl_ipolygon, spl_polygon)

    # Verts

    def get_vert(self, ivert):
        """
        Get ivert vert
        >>> g = Geom((-1,-1,0, 1,-1,0, 0,1,0, 0,0,1), \
                     ((2,1,0), (0,1,3), (1,2,3), (2,0,3)) )  # Good tet
        >>> g.get_vert(2)
        Vector(0.000, 1.000, 0.000)
        """
        return Vector(self.verts[3*ivert:3*ivert+3])

    def append_vert(self, vert):
        """
        Append a vert to the Geom, return its index.
        >>> g = Geom((-1,-1,0, 1,-1,0, 0,1,0, 0,0,1), \
                     ((2,1,0), (0,1,3), (1,2,3), (2,0,3)) )  # Good tet
        >>> g.append_vert((0,0,0))
        4
        """
        self.verts.extend(list(vert))
        return self.get_nverts()-1

    def update_vert(self, ivert, vert):
        """
        Update a vert of a Geom, return its index.
        >>> g = Geom((-1,-1,0, 1,-1,0, 0,1,0, 0,0,1), \
                     ((2,1,0), (0,1,3), (1,2,3), (2,0,3)) )  # Good tet
        >>> g.update_vert(2,(0,0,0)); g.get_vert(2)
        2
        Vector(0.000, 0.000, 0.000)
        """
        self.verts[3*ivert:3*ivert+3] = array.array('f', vert)
        return ivert

    def get_nverts(self):
        """
        Get the len of vertices
        >>> g = Geom((-1,-1,0, 1,-1,0, 0,1,0, 0,0,1), \
                     ((2,1,0), (0,1,3), (1,2,3), (2,0,3)) )  # Good tet
        >>> g.get_nverts()
        4
        """
        return int(len(self.verts)/3)

    def get_iverts(self):
        """
        Get the range for verts
        >>> g = Geom((-1,-1,0, 1,-1,0, 0,1,0, 0,0,1), \
                     ((2,1,0), (0,1,3), (1,2,3), (2,0,3)) )  # Good tet
        >>> g.get_iverts()
        [0, 1, 2, 3]
        """
        return [i for i in range(int(len(self.verts)/3))]

    def merge_duplicated_verts(self, multiplier=1.):
        """
        Remove dup verts, and relink all polygons. No mod of polygons.
        >>> g = Geom((-1,-1,0, 1,-1,0, 0,1,0, 0,0,1, 0,1,0, 0,1,0, \
                      1,-1,0, 1,-1,0,), \
                     ((2,6,0), (0,1,3), (7,4,3), (5,0,3)) )  # Dup verts, 8>4
        >>> g.merge_duplicated_verts(); g  # doctest: +NORMALIZE_WHITESPACE
        Dup verts removed: 4
        4
        Geom(
            (-1.000,-1.000,0.000,  1.000,-1.000,0.000,  0.000,1.000,0.000,
             0.000,0.000,1.000),
            [[2, 1, 0], [0, 1, 3], [1, 2, 3], [2, 0, 3]],
            )
        """
        # Find unique verts and build ivert_to_ivert dict
        original_nverts = self.get_nverts()
        unique_pyverts = []
        ivert_to_ivert = {}
        nverts = self.get_nverts()
        for ivert in range(nverts):
            vert = self.get_vert(ivert)
            seen = False
            for i, selected_pyvert in enumerate(unique_pyverts):
                if (selected_pyvert - vert).is_zero(multiplier):
                    seen = True
                    ivert_to_ivert[ivert] = i
                    break
            if not seen:
                unique_pyverts.append(vert)
                ivert_to_ivert[ivert] = len(unique_pyverts) - 1
        # Update face ivert links
        for i, polygon in enumerate(self.polygons):
            self.update_polygon(
                    i, [ivert_to_ivert[ivert] for ivert in polygon]
                    )
        # Flatten unique_pyverts to build new self.verts
        verts = array.array('f')
        for pyvert in unique_pyverts:
            verts.extend(pyvert)
        self.verts = verts
        print("Dup verts removed:", original_nverts - self.get_nverts())
        return original_nverts - self.get_nverts()

    # Edges

    def _collapse_edge(self, ivert0, ivert1, halfedges):
        # Replace ivert0 with ivert1 in all polygons
        broken_ipolygons = []
        for ipolygon, polygon in enumerate(self.polygons):
            if ivert1 in polygon:
                # Polygon has ivert1 already, try to remove ivert0 only
                try:
                    polygon.remove(ivert0)
                except ValueError:
                    pass
                else:
                    if len(polygon) < 3:
                        broken_ipolygons.append(ipolygon)
            else:
                # Polygon does not have ivert1, try to replace ivert0
                try:
                    polygon[polygon.index(ivert0)] = ivert1
                except ValueError:
                    pass
        # Update halfedges
        for halfedge in halfedges:
            try:
                halfedge[halfedge.index(ivert0)] = ivert1
            except ValueError:
                pass
        # Remove broken ipolygons
        for ipolygon in broken_ipolygons:
            self.remove_polygon(ipolygon)

    def collapse_short_edges(self, limit):  # FIXME develop
        halfedges = self.get_double_halfedges()
        counter = 0
        # Get short halfedges
        limit = limit ** 2
        short_halfedges = []
        for halfedge in halfedges:
            ivert0, ivert1 = halfedge
            v0 = self.get_vert(ivert0)
            v1 = self.get_vert(ivert1)
            if v1.minus(v0).squarelength() < limit:
                short_halfedges.append(list(halfedge))
        # Collapse them
        while short_halfedges:
            ivert0, ivert1 = short_halfedges.pop()
            self._collapse_edge(ivert0, ivert1, short_halfedges)
            counter += 1
        print("Short edges removed:", counter)

    def _split_edge(self, length):  # FIXME develop
        pass

    def split_long_edges(self, length):  # FIXME develop
        pass

    def get_halfedges(self, ipolygons=None):
        """
        Get halfedges dict of ipolygons subset.
        halfedges are: {(1,2):7]} with {(ivert0, ivert1): ipolygon on the left}
        according to iface0 normal up
        >>> g = Geom((-1,-1,0, 1,-1,0, 0,1,0, 0,0,1), \
                     ((2,1,0), (0,1,3), (1,2,3), (2,0,3)) )  # Good tet
        >>> g.get_halfedges(ipolygons=None) # doctest: +NORMALIZE_WHITESPACE
        {(2, 1): 0, (1, 0): 0, (0, 2): 0, (0, 1): 1, (1, 3): 1, (3, 0): 1,
         (1, 2): 2, (2, 3): 2, (3, 1): 2, (2, 0): 3, (0, 3): 3, (3, 2): 3}
        >>> g.get_halfedges(ipolygons=(1,2,3)) # doctest: +NORMALIZE_WHITESPACE
        {(0, 1): 1, (1, 3): 1, (3, 0): 1, (1, 2): 2, (2, 3): 2, (3, 1): 2,
         (2, 0): 3, (0, 3): 3, (3, 2): 3}
        >>> g = Geom((-1,-1,0, 1,-1,0, 0,1,0, 0,0,1), \
                     ((2,1,0), (0,1,3), (1,2,3), (3,0,2)) )  # Unorient tet
        >>> g.get_halfedges()
        Traceback (most recent call last):
        ...
        Exception: ('Non-manifold or unorientable at ipolygon:', 3)
        """
        if ipolygons is None:
            ipolygons = self.get_ipolygons()
        halfedges = dict()
        for ipolygon in ipolygons:
            polygon = self.get_polygon(ipolygon)
            polygon_nverts = len(polygon)
            for i in range(polygon_nverts):
                halfedge = (polygon[i], polygon[(i+1) % polygon_nverts])
                if halfedge in halfedges:
                    raise Exception(
                            'Non-manifold or unorientable at ipolygon:',
                            ipolygon)
                halfedges[halfedge] = ipolygon
        return halfedges

    def get_double_halfedges(self, ipolygons=None):
        """
        Get double halfedges dict of ipolygons subset.
        double halfedges are: {(1,2):7,6]}
        with {(ivert0, ivert1): ipolygon_sx, ipolygon_dx}
        according to iface0 normal up
        >>> g = Geom((-1,-1,0, 1,-1,0, 0,1,0, 0,0,1), \
                     ((2,1,0), (0,1,3), (1,2,3), (2,0,3)) )  # Good tet
        >>> g.get_double_halfedges() # doctest: +NORMALIZE_WHITESPACE
        {(3, 2): (3, 2), (0, 3): (3, 1), (2, 0): (3, 0), (3, 1): (2, 1),
         (1, 2): (2, 0), (0, 1): (1, 0)}
        """
        halfedges = self.get_halfedges(ipolygons)
        halfedges_double = dict()
        while halfedges:
            halfedge, ipolygon_sx = halfedges.popitem()
            opposite = halfedge[1], halfedge[0]
            ipolygon_dx = halfedges.get(opposite, None)
            if ipolygon_dx is not None:
                del halfedges[opposite]
            halfedges_double[halfedge] = ipolygon_sx, ipolygon_dx
        return halfedges_double

    def get_border_halfedges(self, ipolygons=None):
        """
        Get border halfedges dict
        Eg: {(1,2):7]} with {(ivert0, ivert1): ipolygon on the left}
        according to iface0 normal up
        >>> g = Geom((-1,-1,0, 1,-1,0, 0,1,0, 0,0,1), \
                     ((0,1,3), (1,2,3), (2,0,3)) )  # Open tet
        >>> g.get_border_halfedges()
        {(2, 0): 2, (1, 2): 1, (0, 1): 0}
        """
        halfedges = self.get_halfedges(ipolygons)
        border_halfedges = {}
        while halfedges:
            halfedge, ipolygon = halfedges.popitem()
            opposite = halfedge[1], halfedge[0]
            if opposite in halfedges:
                del halfedges[opposite]
            else:
                border_halfedges[halfedge] = ipolygon
        return border_halfedges

    def get_multi_halfedges(self, ipolygons=None):
        """
        Get multi halfedges, that are common between polygons
        Eg: {(1,2,3,4):[7,4]} with
        {(ivert0, ivert1): [ipolygon on the left, ipolygon on the right]}
        >>> g = Geom((0,0,0, 1,0,0, 2,0,0, 3,0,0, 4,0,0, 5,0,0, 6,0,0,\
                      3,1,0, 3,3,0, 3,-3,0),\
                     ((0,1,2,3,7,4,5,6,8), (9,6,5,4,3,2,1,0), (4,7,3)) )
        >>> g.get_multi_halfedges()  # doctest: +NORMALIZE_WHITESPACE
        {(6, 8, 0): (0, None), (6, 5, 4): (1, 0), (3, 2, 1, 0): (1, 0),
         (0, 9, 6): (1, None), (4, 7, 3): (2, 0)}
        """
        halfedges = self.get_halfedges(ipolygons)
        # Get polygons common halfedges
        # {(ipolygon_sx, ipolygon_dx): [halfedge, halfedge, ...]
        common_halfedges = {}
        while halfedges:
            halfedge, ipolygon_sx = halfedges.popitem()
            opposite = halfedge[1], halfedge[0]
            if opposite in halfedges:
                ipolygon_dx = halfedges[opposite]
                del halfedges[opposite]
                try:
                    common_halfedges[
                            (ipolygon_sx, ipolygon_dx)
                            ].append(halfedge)
                except KeyError:
                    common_halfedges[(ipolygon_sx, ipolygon_dx)] = [halfedge, ]
            else:  # it is a border
                try:
                    common_halfedges[(ipolygon_sx, None)].append(halfedge)
                except KeyError:
                    common_halfedges[(ipolygon_sx, None)] = [halfedge, ]
        # Get polygons multi halfedges
        multi_halfedges = {}
        while common_halfedges:
            ipolygons, halfedges = common_halfedges.popitem()
            while halfedges:
                multi = list(halfedges.pop())
                done = False
                while not done:
                    done = True
                    for candidate in halfedges:
                        if candidate[0] == multi[-1]:
                            multi.append(candidate[1])
                            done = False
                            halfedges.remove(candidate)
                            break
                        elif candidate[1] == multi[0]:
                            multi.insert(0, candidate[0])
                            done = False
                            halfedges.remove(candidate)
                            break
                if len(multi) > 2:
                    multi_halfedges[tuple(multi)] = ipolygons
        return multi_halfedges

    def remove_multi_halfedges(self):  # FIXME Test
        """
        Remove multi halfedges from polygon borders, retaining manifoldness
        """
        multi_halfedges = self.get_multi_halfedges()
        for multi, ipolygons in multi_halfedges.items():
            polygon0 = self.get_polygon(ipolygons[0])
            if ipolygons[1] is not None:  # An edge between coplanar polygons
                polygon1 = self.get_polygon(ipolygons[1])
                for m in multi[1:-1]:
                    polygon0.remove(m)
                    polygon1.remove(m)
                self.update_polygon(ipolygons[0], polygon0)
                self.update_polygon(ipolygons[1], polygon1)
            else:  # An external edge, remove only if collinear
                nverts = len(multi)-2
                for i in range(nverts):
                    v0 = self.get_vert(multi[i])
                    p = self.get_vert(multi[i+1])
                    v1 = self.get_vert(multi[i+2])
                    if p.is_within(v0, v1) and p.is_collinear(
                            v0, v1,
                            multiplier=2
                            ):
                        polygon0.remove(multi[i+1])
                self.update_polygon(ipolygons[0], polygon0)

    def get_border_loops(self, ipolygons):
        """
        Get oriented border vert loops,
        Eg: {(ivert0, ivert1, ...):[ipolygon, ...]}
        according to iface0 normal up
        >>> g = Geom((-1,-1,0, 1,-1,0, 1,1,0, -1,1,0, -3, 1,0, -3,-1,0, \
               3,-1,0, 3, 1,0, 1,3,0, -1,3,0, -1,-3,0,  1,-3,0),\
             ((0,1,2,3), (5,0,3,4), (1,6,7,2), (3,2,8,9), (10,11,1,0)),\
             (0,1,2,3,4), \
            )  # Open clover on z=0, n=+k
        >>> g.get_border_loops(ipolygons=(0,1,2,3,4,))
        {(5, 0, 10, 11, 1, 6, 7, 2, 8, 9, 3, 4): [0, 1, 2, 3, 4]}
        """
        border_halfedges = self.get_border_halfedges(ipolygons)
        # Get loops
        loops = {}
        while border_halfedges:
            halfedge, ipolygon = border_halfedges.popitem()
            loop = [halfedge[0], halfedge[1]]
            loop_ipolygons = [ipolygon, ]
            # While loop not closed
            while loop[-1] != loop[0]:
                # Get candidates for the loop
                candidates = []
                for halfedge, ipolygon in border_halfedges.items():
                    if loop[-1] == halfedge[0]:
                        candidates.append((halfedge, ipolygon))
                # Good vert or singularity?
                # Choose the candidate that goes home
                if len(candidates) == 1:
                    halfedge, ipolygon = candidates[0]
                elif len(candidates) == 2:
                    ve = Vector(self.get_vert(loop[-2])) \
                         - Vector(self.get_vert(loop[-1]))
                    v0 = Vector(self.get_vert(candidates[0][0][0])) \
                        - Vector(self.get_vert(candidates[0][0][1]))
                    v1 = Vector(self.get_vert(candidates[1][0][0])) \
                        - Vector(self.get_vert(candidates[1][0][1]))
                    dot0 = v0.dot(ve)
                    dot1 = v1.dot(ve)
                    if dot0 < dot1:
                        halfedge, ipolygon = candidates[0]
                    else:
                        halfedge, ipolygon = candidates[1]
                else:  # In case of 0 or more than 2 candidates
                    raise Exception('Invalid GEOM, non-manifold in:',
                                    self.get_vert(loop[-1]))
                # Append candidate
                del(border_halfedges[halfedge])
                loop.append(halfedge[1])
                loop_ipolygons.append(ipolygon)
            # Closed
            loop.pop()  # Remove closing ivert that is duplicated
            loops[tuple(loop)] = loop_ipolygons
            # Get all included polygons
            for loop, loop_ipolygons in loops.items():
                loops[loop] = self.get_cont_ipolygons(
                        loop_ipolygons[0], ipolygons
                        )
        return loops

    def get_bordering_ipolygons(self, ipolygon, ipolygons=None):
        """
        Get a set of polygons in ipolygons that border ipolygon
        >>> g = Geom((-1,-1,0, 1,-1,0, 1,1,0, -1,1,0, -3, 1,0, -3,-1,0, \
               3,-1,0, 3, 1,0, 1,3,0, -1,3,0, -1,-3,0,  1,-3,0),\
             ((0,1,2,3), (5,0,3,4), (1,6,7,2), (3,2,8,9),\
             (10,11,1,0), (4,3,9), (8,2,7), (6,1,11), (0,5,10),) \
            )  # Open clover w addition on z=0, n=+k
        >>> #  y ↑
        >>> #   9─8
        >>> #  /│ │ \
        >>> # 4─3─2─7
        >>> # │ │∙│ │→ x
        >>> # 5─0─1─6
        >>> #  \│ │/
        >>> #  10-11
        >>> g.get_bordering_ipolygons(ipolygon=0, ipolygons=None)
        {1, 2, 3, 4}
        """
        bordering_ipolygons = []
        ipolygons_halfedges = self.get_halfedges(ipolygons)
        for halfedge in self.get_halfedges((ipolygon,)):
            try:
                bordering_ipolygons.append(
                        ipolygons_halfedges[(halfedge[1], halfedge[0])]
                        )
            except KeyError:
                pass
        return set(bordering_ipolygons)

    def get_cont_ipolygons(self, ipolygon, ipolygons):
        """
        Get a list of polygons in ipolygons that are continuous to ipolygon
        >>> g = Geom((-1,-1,0, 1,-1,0, 1,1,0, -1,1,0, -3, 1,0, -3,-1,0, \
               3,-1,0, 3, 1,0, 1,3,0, -1,3,0, -1,-3,0,  1,-3,0),\
             ((0,1,2,3), (5,0,3,4), (1,6,7,2), (3,2,8,9),\
             (10,11,1,0), (4,3,9), (8,2,7), (6,1,11), (0,5,10),) \
            )  # Open clover w addition on z=0, n=+k
        >>> #  y ↑
        >>> #   9─8
        >>> #  /│ │ \
        >>> # 4─3─2─7
        >>> # │ │∙│ │→ x
        >>> # 5─0─1─6
        >>> #  \│ │/
        >>> #  10-11
        >>> g.get_cont_ipolygons(ipolygon=0, ipolygons=(0,1,2,3,4,5,6,7,8,))
        [0, 1, 2, 3, 4, 5, 6, 7, 8]
        """
        cont_ipolygons = set([ipolygon, ])
        new_ipolygons = set([ipolygon, ])
        ipolygons = set(ipolygons)
        while new_ipolygons:
            ipolygon = new_ipolygons.pop()
            bord = self.get_bordering_ipolygons(ipolygon, ipolygons)
            new_ipolygons |= bord
            ipolygons -= bord
            cont_ipolygons |= bord
        return list(cont_ipolygons)

    def _get_tri_earclip(self, polygon, normal):
        """
        Get valid earclip of polygon, remove it from polygon, return it
        """
        polygon_nverts = len(polygon)
        for i0 in range(polygon_nverts):
            i1, i2 = (i0+1) % polygon_nverts, (i0+2) % polygon_nverts
            ivert0, ivert1, ivert2 = polygon[i0], polygon[i1], polygon[i2]
            a, b, c = \
                self.get_vert(ivert0),\
                self.get_vert(ivert1),\
                self.get_vert(ivert2)
            # Test counter-clockwise ear
            ccw = b.minus(a).cross(c.minus(b)).dot(normal) > 0.
            # Test no other vert in the ear
            vwt = any((self.get_vert(p).is_within_tri(a, b, c, normal)
                       for j, p in enumerate(polygon)
                       if j not in (i0, i1, i2)
                       ))
            # If ok, send ear
            if ccw and not vwt:
                del(polygon[i1])
                return polygon, (ivert0, ivert1, ivert2)
        raise Exception('Triangulation impossible, tri:', a, b, c, normal)

    def _get_convex_earclip(self, polygon, normal, tri=False):
        pass # FIXME stub

    def get_tris_of_polygon(self, ipolygon):
        """
        Triangulate ipolygon with no zero-area tris
        >>> g = Geom((0,0,0, 3,0,0, 3,1,0, 1,1,0, 1,3,0, 0,3,0,), \
                     ((5,0,1,2,3,4,), ))    # L concave
        >>> g.get_tris_of_polygon(ipolygon=0)  # doctest: +NORMALIZE_WHITESPACE
        [(0, 1, 2), (0, 2, 3), (5, 0, 3), (5, 3, 4)]
        >>> g = Geom((0,0,0, 1,0,0, 2,0,0, 3,0,0, 3,1,0, 2,1,0, 1,1,0, 1,2,0,\
                      1,3,0, 0,3,0, 0,2,0, 0,1,0), \
                     ((0,1,2,3,4,5,6,7,8,9,10,11), ))  # L concave, alignment
        >>> g.get_tris_of_polygon(ipolygon=0)  # doctest: +NORMALIZE_WHITESPACE
        [(2, 3, 4), (1, 2, 4), (0, 1, 4), (0, 4, 5), (0, 5, 6), (0, 6, 7),
         (0, 7, 8), (8, 9, 10), (8, 10, 11), (0, 8, 11)]
        >>> g = Geom((0,0,0, 1,0,0, 2,0,0, 2,1,0, ), \
                     ((0,1,2,3,2), ))    # Zero area
        >>> g.get_tris_of_polygon(ipolygon=0)  # doctest: +NORMALIZE_WHITESPACE
        Traceback (most recent call last):
        ...
        Exception: ('Triangulation impossible, tri:',
                    Vector(2.000, 0.000, 0.000), Vector(0.000, 0.000, 0.000),
                    Vector(1.000, 0.000, 0.000), Vector(0.000, 0.000, -1.000))
        >>> g = Geom((0,0,0, 1,0,0, 1,0,0, 1,2,0, ), \
                     ((0,1,2,3,), ))    # Zero edge
        >>> g.get_tris_of_polygon(ipolygon=0)  # doctest: +NORMALIZE_WHITESPACE
        Traceback (most recent call last):
        ...
        Exception: ('Triangulation impossible, tri:',
                    Vector(1.000, 2.000, 0.000), Vector(0.000, 0.000, 0.000),
                    Vector(1.000, 0.000, 0.000), Vector(0.000, 0.000, 1.000))
        """
        polygon = self.get_polygon(ipolygon)[:]
        polygon_nverts = len(polygon)
        # Short cut
        if polygon_nverts == 3:
            return [tuple(polygon), ]
        # Get the polygon overall normal
        normal = self.get_polygon_normal(ipolygon)
        # Search for triangulation
        tris = []
        while len(polygon) > 2:
            polygon, tri = self._get_tri_earclip(polygon, normal)
            tris.append(tri)
        return tris

    # STL/OBJ

    def to_STL(self, filepath):
        """
        Write self to STL file
        >>> g = Geom((-1.0, -1.0, -1.0,  -1.0, -1.0, 1.0,  -1.0, 1.0,  1.0, \
                      -1.0,  1.0, -1.0,   1.0,  1.0, 1.0,   1.0, 1.0, -1.0, \
                       1.0, -1.0, -1.0,   1.0, -1.0, 1.0), \
                    ((0,1,2,3), (7,6,5,4), (1,0,6,7), (2,4,5,3), \
                     (0,3,5,6), (1,7,4,2)),  \
                     ('Red', 'Magenta', 'Green', 'Yellow', 'Blue', 'Cyan',)) \
                    # A good cube w surfid
        >>> g.to_STL('../test/doctest.stl')
        to_STL: ../test/doctest.stl
        >>> Geom.from_STL('../test/doctest.stl') \
            # doctest: +NORMALIZE_WHITESPACE
        Dup verts removed: 28
        Geom(
            (-1.000,-1.000,-1.000,  -1.000,-1.000,1.000,  -1.000,1.000,1.000,
             -1.000,1.000,-1.000,  1.000,-1.000,1.000,  1.000,-1.000,-1.000,
             1.000,1.000,-1.000,  1.000,1.000,1.000),
            [[0, 1, 2], [0, 2, 3], [4, 5, 6], [4, 6, 7], [1, 0, 5], [1, 5, 4],
             [2, 7, 6], [2, 6, 3], [0, 3, 6], [0, 6, 5], [1, 4, 7], [1, 7, 2]],
            )
        >>> g.to_STL('../test/doctest.stl')
        to_STL: ../test/doctest.stl
        """
        with open(filepath, 'w') as f:
            f.write('solid name\n')
            for ipolygon in range(self.get_npolygons()):
                tris = self.get_tris_of_polygon(ipolygon)
                for tri in tris:
                    f.write('facet normal 0 0 0\n')
                    f.write(' outer loop\n')
                    for ivert in tri:
                        f.write('  vertex '
                                '{v[0]:.9f} {v[1]:.9f} {v[2]:.9f}\n'.format(
                                  v=self.get_vert(ivert)),
                                )
                    f.write(' endloop\n')
                    f.write('endfacet\n')
            f.write('endsolid name\n')
        print('to_STL:', filepath)

    @classmethod
    def from_STL(cls, filepath, surfid='White'):
        """
        Get new Geom from STL file
        Doctest in Geom.to_STL()
        """
        import os
        path, filename = os.path.split(filepath)
        # Get STL mesh
        from stl import mesh
        mesh = mesh.Mesh.from_file(filepath)
        verts, polygons, surfids, py_verts = [], [], [], []
        for iface, p in enumerate(mesh.points):
            # p is [-1.0, -1.0, -1.0, -1.0, -1.0, 1.0, -1.0, 1.0, 1.0]
            verts.extend(p)
            polygons.append((3*iface, 3*iface+1, 3*iface+2))
            py_verts.append((p[0], p[1], p[2]))
            py_verts.append((p[3], p[4], p[5]))
            py_verts.append((p[6], p[7], p[8]))
            surfids.append(surfid)
        g = Geom(verts, polygons, surfids, hid=filename)
        g.merge_duplicated_verts()
        g.check_geom_sanity()
        return g

    def to_OBJ(self, filepath, triangulate=False):
        """
        Write self to OBJ file
        >>> g = Geom((-1.0, -1.0, -1.0,  -1.0, -1.0, 1.0,  -1.0, 1.0,  1.0, \
                      -1.0,  1.0, -1.0,   1.0,  1.0, 1.0,   1.0, 1.0, -1.0, \
                       1.0, -1.0, -1.0,   1.0, -1.0, 1.0), \
                    ((0,1,2,3), (7,6,5,4), (1,0,6,7), (2,4,5,3), \
                     (0,3,5,6), (1,7,4,2)),  \
                     ('Red', 'Magenta', 'Green', 'Yellow', 'Blue', 'Cyan',)) \
                    # A good cube w surfid
        >>> g.to_OBJ('../test/doctest.obj')
        to_OBJ: ../test/doctest.obj
        >>> g = Geom.from_OBJ('../test/doctest.obj').popitem()[1]
        from_OBJ: ../test/doctest.obj -> doctest.obj
        Dup verts removed: 0
        >>> g.get_polygon(1); g.get_vert(7)
        [7, 6, 5, 4]
        Vector(1.000, -1.000, 1.000)
        >>> g.to_OBJ('../test/doctest2.obj', triangulate=True)
        to_OBJ: ../test/doctest2.obj
        """
        import os
        path, filename = os.path.split(filepath)
        # Arrange polygons by surfid
        surfid_to_ipolygons = {}  # FIXME move to def, duplicated
        for ipolygon, _ in enumerate(self.polygons):
            surfid = self.get_polygon_surfid(ipolygon)
            try:
                surfid_to_ipolygons[surfid].append(ipolygon)
            except KeyError:
                surfid_to_ipolygons[surfid] = [ipolygon, ]
        # Write geometry
        with open(filepath, 'w') as f:
            f.write('mtllib default.mtl\n')
            f.write('o {}\n'.format(self.hid or filename))
            for ivert in self.get_iverts():
                vert = self.get_vert(ivert)
                new_vert = (vert[0], vert[2], -vert[1])  # Different ref sys
                f.write('v {0[0]} {0[1]} {0[2]}\n'.format(new_vert))
            for surfid, ipolygons in surfid_to_ipolygons.items():
                f.write('usemtl {}\n'.format(surfid))
                for ipolygon in ipolygons:
                    if triangulate:
                        polygons = self.get_tris_of_polygon(ipolygon)
                    else:
                        polygons = (self.get_polygon(ipolygon),)
                    for polygon in polygons:
                        str_polygon = ' '.join(
                                [str(ivert+1) for ivert in polygon]
                                )
                        f.write('f {}\n'.format(str_polygon))
        # Write predefined materials
        with open('{}/default.mtl'.format(path), 'w') as f:
            f.write(
                    """
                    # Materials
                    newmtl White\nKd 0.6 0.6 0.6
                    newmtl Red\nKd 0.6 0.0 0.0
                    newmtl Green\nKd 0.0 0.6 0.0
                    newmtl Blue\nKd 0.0 0.0 0.6
                    newmtl Yellow\nKd 0.6 0.6 0.0
                    newmtl Cyan\nKd 0.0 0.6 0.6
                    newmtl Magenta\nKd 0.6 0.0 0.6
                    """
                    )
        print('to_OBJ:', filepath)

    @classmethod
    def from_OBJ(cls, filepath):
        """
        Get a dict of new Geom from OBJ file
        Doctest in Geom.to_OBJ()
        """
        # Init
        import os
        path, filename = os.path.split(filepath)
        geoms = {}
        current_geom = Geom()
        current_surfid = 'White'
        nverts = 0
        # Read file
        with open(filepath, 'r') as f:
            for line in f:
                tokens = line[:-1].split(' ')
                # vert
                if tokens[0] == 'v':
                    current_geom.append_vert((
                            float(tokens[1]),
                            -float(tokens[3]),
                            float(tokens[2]),
                            ))  # Different ref sys
                # polygon
                elif tokens[0] == 'f':
                    # Remove texture and normal info
                    tokens = [t.split('/')[0] for t in tokens]
                    # Check polygon
                    polygon = [int(t)-1-nverts for t in tokens[1:]]
                    for ivert in polygon:
                        if ivert < 0:
                            raise Exception(
                                    'OBJ format not supported, negative ivert.'
                                    )
                    # Get polygon
                    current_geom.append_polygon(
                            polygon=polygon,
                            surfid=current_surfid,
                            normal=(1, 0, 0),  # updated later
                            )
                # surfid
                elif tokens[0] == 'usemtl':
                    current_surfid = tokens[1]
                # geom
                elif tokens[0] == 'o':
                    if not geoms:
                        geoms[tokens[1]] = current_geom
                        current_geom.hid = tokens[1]
                    else:
                        nverts = current_geom.get_nverts()
                        current_geom = Geom(hid=tokens[1])
                        geoms[tokens[1]] = current_geom
                        current_surfid = 'White'

        if not geoms:
            geoms[filename] = current_geom
        for hid, g in geoms.items():
            print('from_OBJ:', filepath, '->', hid)
            g.update_normals()
            g.merge_duplicated_verts()
            g.check_geom_sanity()
        return geoms

    # Geometry sanity

    def check_loose_verts(self):
        """
        Check loose vertices: vertices that have no connectivity
        >>> g = Geom((-1,-1,0, 1,-1,0, 0,1,0, 0,0,1, 0,0,0), \
                     ((2,1,0), (0,1,3), (1,2,3), (2,0,3)) )  # Loose vert
        >>> g.check_loose_verts()
        Traceback (most recent call last):
        ...
        Exception: Invalid GEOM, loose verts.
        """
        nverts = self.get_nverts()
        npolygons = self.get_npolygons()
        used_iverts = []
        for ipolygon in range(npolygons):
            used_iverts.extend(self.get_polygon(ipolygon))
        used_iverts = set(used_iverts)
        if nverts != len(used_iverts) or nverts != max(used_iverts) + 1:
            raise Exception("Invalid GEOM, loose verts.")

    def check_degenerate_geometry(self):
        """
        Check degenerate geometry, as zero lenght edges and zero area faces
        """
        for ipolygon in range(self.get_npolygons()):
            self.get_tris_of_polygon(ipolygon)

    def check_is_solid(self):
        """
        Check surface:
        - 2-manifold and closed, each edge should join two faces,
        - orientable, adjoining faces should have same normals
        >>> g = Geom((-1,-1,0, 1,-1,0, 0,1,0, 0,0,1), \
                     ((0,1,3), (1,2,3), (2,0,3)) )  # Open tet
        >>> g.check_is_solid()
        Traceback (most recent call last):
        ...
        Exception: ('Non closed at ipolygons:', [2, 1, 0])
        """
        border_halfedges = self.get_border_halfedges()
        if border_halfedges:
            raise Exception("Non closed at ipolygons:",
                            [b for b in border_halfedges.values()])

#def check_euler(self):  # FIXME test
#    """
#    Euler formula: nverts - nedges + nfaces = chi
#    Euler characteristic chi of the connected sum of g tori is:
#    chi = 2 − 2g, with g genus
#    g = 0, 1, 2, 3, ... => chi = 2, 0, -2, -4, ...
#    """
#    nverts = get_nverts(igeom)
#    nfaces = get_nfaces(igeom)
#    # Count edges
#    nedges = 0
#    seen_halfedges = []
#    halfedges = get_halfedges(igeom, get_ifaces(igeom))
#    for halfedge in halfedges:
#        opposite = halfedge[1], halfedge[0]
#        if halfedge in seen_halfedges or opposite in seen_halfedges:
#            continue
#        else:
#            seen_halfedges.append(halfedge)
#            nedges += 1
#    # Check Euler formula
#    chi = nverts - nedges + nfaces
#    if chi not in range(2, 100, 2):
#        raise Exception('Invalid GEOM, chi in Euler formula is:', chi)

    def check_geom_sanity(self):
        """
        Check geometry sanity

        If the mesh is correct and encloses a volume, this can be checked with
        prior tests: checking orientability, non-borders,
        non-self-intersecting.
        After that we can calculate its topological features and check if
        Euler's formula  C+V=A+2(S-H) is satisfied.
        If the mesh is not correct, many geometric algorithms will fail.
        The only solution in this case is the user repairing the mesh.
        >>> g = Geom((-1.0, -1.0, -1.0,  -1.0, -1.0, 1.0,  -1.0, 1.0,  1.0, \
                      -1.0,  1.0, -1.0,   1.0,  1.0, 1.0,   1.0, 1.0, -1.0, \
                       1.0, -1.0, -1.0,   1.0, -1.0, 1.0), \
                    ((0,1,2,3), (7,6,5,4), (1,0,6,7), (2,4,5,3), \
                     (0,3,5,6), (1,7,4,2)),  \
                     ('Red', 'Magenta', 'Green', 'Yellow', 'Blue', 'Cyan',)) \
                    # A good cube w surfid
        >>> g.check_geom_sanity()
        """
        self.check_loose_verts()
        self.check_degenerate_geometry()
        self.check_is_solid()
#        self.check_euler()
#        self.check_flat_polygons()
        # Check correct normals for a solid in fluid FIXME working here
        # Check self intersection  # FIXME not working

    # Boolean

    def union(self, geom):  # FIXME
        """
        Update current geom to represent union with geom.
        """
        pass


# BSP tree

class BSPNode(object):
    """
    A node in a BSP tree.
    >>> g = Geom((-1,-1,0, 1,-1,0, 0,1,0, 0,0,1), \
                 ((2,1,0), (0,1,3), (1,2,3), (2,0,3)) )  # Good tet
    >>> n = BSPNode(geom=g); n.build(); n  # doctest: +NORMALIZE_WHITESPACE
    BSP tree - geom hid: None, ipolygons: [0]
        │ plane: Plane(normal=Vector(0.000, 0.000, -1.000), distance=0.000)
        ├─front_node: None
        └─back_node: geom hid: None, ipolygons: [1]
          │ plane: Plane(normal=Vector(0.000, -0.707, 0.707), distance=0.707)
          ├─front_node: None
          └─back_node: geom hid: None, ipolygons: [2]
            │ plane: Plane(normal=Vector(0.816, 0.408, 0.408), distance=0.408)
            ├─front_node: None
            └─back_node: geom hid: None, ipolygons: [3]
              │ plane: Plane(normal=Vector(-0.816, 0.408, 0.408),
                             distance=0.408)
              ├─front_node: None
              └─back_node: None
    >>> g = Geom((0,1,0, 1,1,0, 1,2,0, -1,2,0, -1,-2,0, 1,-2,0, 1,-1,0, \
                  0,-1,0,\
                  0,1,1, 1,1,1, 1,2,1, -1,2,1, -1,-2,1, 1,-2,1, 1,-1,1, \
                  0,-1,1,\
                  ), \
                 ((3,2,1,0), (7,4,3,0), (7,6,5,4), \
                  (8,9,10,11), (8,11,12,15), (12,13,14,15), \
                  (12,11,3,4), (7,0,8,15), (9,8,0,1), (6,7,15,14), \
                  (11,10,2,3), (13,12,4,5), (1,2,10,9), (5,6,14,13), \
                  ) )  # Good concave C shape
    >>> n = BSPNode(geom=g); n.build(); n  # doctest: +NORMALIZE_WHITESPACE
    BSP tree - geom hid: None, ipolygons: [0, 1, 2]
        │ plane: Plane(normal=Vector(0.000, 0.000, -1.000), distance=0.000)
        ├─front_node: None
        └─back_node: geom hid: None, ipolygons: [3, 4, 5]
          │ plane: Plane(normal=Vector(0.000, 0.000, 1.000), distance=1.000)
          ├─front_node: None
          └─back_node: geom hid: None, ipolygons: [6]
            │ plane: Plane(normal=Vector(-1.000, 0.000, 0.000), distance=1.000)
            ├─front_node: None
            └─back_node: geom hid: None, ipolygons: [7]
              │ plane: Plane(normal=Vector(1.000, 0.000, 0.000),
                             distance=0.000)
              ├─front_node: geom hid: None, ipolygons: [8]
              │ │ plane: Plane(normal=Vector(0.000, -1.000, 0.000),
                               distance=-1.000)
              │ ├─front_node: geom hid: None, ipolygons: [9]
                │ │ plane: Plane(normal=Vector(0.000, 1.000, 0.000),
                                 distance=-1.000)
                │ ├─front_node: None
                │ └─back_node: geom hid: None, ipolygons: [11]
                    │ plane: Plane(normal=Vector(0.000, -1.000, 0.000),
                                   distance=2.000)
                    ├─front_node: None
                    └─back_node: geom hid: None, ipolygons: [13]
                      │ plane: Plane(normal=Vector(1.000, 0.000, 0.000),
                                     distance=1.000)
                      ├─front_node: None
                      └─back_node: None
              │ └─back_node: geom hid: None, ipolygons: [10]
                  │ plane: Plane(normal=Vector(0.000, 1.000, 0.000),
                                 distance=2.000)
                  ├─front_node: None
                  └─back_node: geom hid: None, ipolygons: [12]
                    │ plane: Plane(normal=Vector(1.000, 0.000, 0.000),
                                   distance=1.000)
                    ├─front_node: None
                    └─back_node: None
              └─back_node: geom hid: None, ipolygons: [14]
                │ plane: Plane(normal=Vector(0.000, 1.000, 0.000),
                               distance=2.000)
                ├─front_node: None
                └─back_node: geom hid: None, ipolygons: [15]
                  │ plane: Plane(normal=Vector(0.000, -1.000, 0.000),
                                 distance=2.000)
                  ├─front_node: None
                  └─back_node: None
    >>> g.to_OBJ('../test/c-shape.obj')
    to_OBJ: ../test/c-shape.obj
    """
    def __init__(self, geom):
        # Tree
        self.plane = None       # Cutting Plane instance
        self.front_node = None  # Front BSPNode, for tree
        self.back_node = None   # Back BSPNode, for tree
        # Geom
        self.geom = geom        # Link to related geom
        self.ipolygons = []     # Coplanar ipolygons

    def __repr__(self):
        return 'BSP tree - {}'.format(self._repr_tree())

    def _repr_tree(self, back=True):
        # Get children trees
        front_tree = "None"
        if self.front_node:
            front_tree = self.front_node._repr_tree(back=False)
        back_tree = "None"
        if self.back_node:
            back_tree = self.back_node._repr_tree(back=True)
        # Join texts
        line = "│ "
        if back:
            line = "  "
        header = 'geom hid: {0}, ipolygons: {1}\n'.format(
            self.geom.hid,
            self.ipolygons,
            )
        text = '{3}│ plane: {0}\n{3}├─front_node: {1}\n'\
            '{3}└─back_node: {2}'.format(
                    self.plane,
                    front_tree,
                    back_tree,
                    line,
                    )
        return header + textwrap.indent(text, '  ')

    def clone(self):
        node = BSPNode(self.geom)
        if self.plane:
            node.plane = self.plane.clone()
        if self.front_node:
            node.front_node = self.front_node.clone()
        if self.back_node:
            node.back = self.back_node.clone()
        node.ipolygons = self.ipolygons[:]
        return node

    def _invert_node(self):
        """
        Invert node orientation, solid <-> void
        """
        # Flip normals
        self.plane.flip()
        # Invert tree
        if self.front_node:
            self.front_node._invert_node()
        if self.back_node:
            self.back_node._invert_node()
        # Swap front and back nodes
        temp = self.front_node
        self.front_node = self.back_node
        self.back_node = temp

    def invert(self):
        """
        Swap solid space and empty space.
        """
        self.geom.flip()  # Invert geometry, once for all
        self._invert_node()  # Invert bsp tree

    def clip_polygons(self, geom, ipolygons):
        """
        Recursively remove all geom ipolygons that are inside this BSP tree.
        """
        if not self.plane:
            return

        # Split all geom ipolygons by self.plane, accumulate fragments
        front = []  # front ipolygons
        back = []   # back ipolygons
        for ipolygon in ipolygons:
            geom.split_polygon(ipolygon, self.plane,
                               front, back, front, back)

        if self.front_node:
            front = self.front_node.clip_polygons(geom, front)

        if self.back_node:
            back = self.back_node.clip_polygons(geom, back)
        else:
            back = []  # Remove back polygons of a leaf without back_node

        front.extend(back)
        return front

    def clip_to(self, clipping_bsp):
        """
        Remove all polygons in this BSP tree that are inside the
        clipping BSP tree
        """
        self.ipolygons = clipping_bsp.clip_polygons(self.geom, self.ipolygons)
        if self.front_node:
            self.front_node.clip_to(clipping_bsp)
        if self.back_node:
            self.back_node.clip_to(clipping_bsp)

    def get_all_ipolygons(self):
        """
        Return a list of all ipolygons in this BSP tree.
        """
        ipolygons = self.ipolygons[:]
        if self.front_node:
            ipolygons.extend(self.front_node.get_all_ipolygons())
        if self.back_node:
            ipolygons.extend(self.back_node.get_all_ipolygons())
        return ipolygons

    def build(self, ipolygons=None):
        """
        Recursively build a BSP tree out of the geom polygons.
        When called on an existing tree, the new polygons (from the same geom)
        are filtered down to the bottom of the tree and become new nodes there.
        The splitting plane is set as the first polygon plane.
        """
        # Protect
        if ipolygons is None:
            ipolygons = self.geom.get_ipolygons()
        if not ipolygons:
            return None

        # Set the cutting plane
        i = 0
        if not self.plane:
            i = 1
            self.plane = self.geom.get_polygon_plane(ipolygons[0])
            self.ipolygons.append(ipolygons[0])

        # Split all polygons using self.plane
        front = []  # front ipolygons
        back = []   # back ipolygons
        for ipolygon in ipolygons[i:]:  # If got for plane, start from the 2nd
            # Coplanar front and back polygons go into self.ipolygons
            self.geom.split_polygon(ipolygon, self.plane,
                                    self.ipolygons, self.ipolygons,
                                    front, back)

        # Recursively build the BSP tree
        if len(front) > 0:
            if not self.front_node:
                self.front_node = BSPNode(self.geom)
            self.front_node.build(front)
        if len(back) > 0:
            if not self.back_node:
                self.back_node = BSPNode(self.geom)
            self.back_node.build(back)

    def sync_geom(self):  # FIXME why not a new Geom?
        """
        Sync polygons from self to self.geom
        """
        geom = self.geom
        new_polygons, new_surfids, new_normals = [], [], []
        new_surfids = []
        for ipolygon in self.get_all_ipolygons():
            new_polygons.append(geom.get_polygon(ipolygon))
            new_surfids.append(geom.get_polygon_surfid(ipolygon))
            new_normals.append(geom.get_polygon_normal(ipolygon))
        geom.polygons = new_polygons
        geom.surfids = new_surfids
        geom.normals = new_normals

    def merge_polygons_to_concave(self):  # FIXME move to geom
        """
        Merge coplanar polygons with same surfid to concave polygons
        BSP tree cannot be used any more at the end.
        """
        self.geom.merge_polygons_to_concave(self.ipolygons)
        # Do the same for all the tree
        if self.front_node:
            self.front_node.merge_polygons_to_concave()
        if self.back_node:
            self.back_node.merge_polygons_to_concave()


if __name__ == "__main__":
    import doctest
    doctest.testmod()

    name = "torus"

    geometries = Geom.from_OBJ('../test/{0}/{0}.obj'.format(name))
    g = geometries['{}_a'.format(name)]
    h = geometries['{}_b'.format(name)]

    # Create and build BSP trees
    a = BSPNode(g)
    a.build()

    b = BSPNode(h)
    b.build()

    # Remove each interior
    a.clip_to(b)
    b.clip_to(a)

    # Remove shared coplanars
    b.invert()
    b.clip_to(a)
    b.invert()

    # Merge coplanar polygons with same surfid
#    a.merge_polygons_to_concave()
#    b.merge_polygons_to_concave()

    # Sync
    a.sync_geom()
    b.sync_geom()

#    g.collapse_short_edges(limit=.001)
#    h.collapse_short_edges(limit=.001)

    g.remove_multi_halfedges()
    h.remove_multi_halfedges()

    g.append(h)
#    g.collapse_short_edges(limit=.01)

#    g.to_OBJ('../test/{0}/{0}_union.obj'.format(name), triangulate=True)
    g.to_OBJ('../test/{0}/{0}_union.obj'.format(name), triangulate=False)
