#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 10 17:12:27 2018

@author: egissi
"""

import math
import array

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
        """
        Linear interpolation from self to a.
        >>> Vector(0,0,1).lerp(Vector(10,0,1),.7)
        Vector(7.000, 0.000, 1.000)
        """
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

    def __eq__(self, other):  # FIXME what isclose?
        return \
            math.isclose(self.x, other.x, rel_tol=EPSILON) and \
            math.isclose(self.y, other.y, rel_tol=EPSILON) and \
            math.isclose(self.z, other.z, rel_tol=EPSILON)

    def __hash__(self):
        return hash((self.x, self.y, self.z))

    def is_zero(self):
        return abs(self.x) < EPSILON and \
               abs(self.y) < EPSILON and \
               abs(self.z) < EPSILON

    def is_collinear(self, b, c):
        """
        Return true if self, b, and c all lie on the same line.
        >>> a, b, c = Vector(0,0,0), Vector(1,0,0), Vector(2,0,0)
        >>> a.is_collinear(b, c)
        True
        >>> a, b, c = Vector(0,0,0), Vector(0,-1,0), Vector(0,-2,1)
        >>> a.is_collinear(b, c)
        False
        >>> a, b, c = Vector(0,0,0), Vector(0,0,0), Vector(0,0,0)
        >>> a.is_collinear(b, c)
        True
        """
        # Proposed elsewhere, to correctly use epsilon FIXME
        # area_tri = abs((c-self).cross(b-self).length())
        # area_square = max((c-self).length() ** 2,
        #   (b-self).length() ** 2, (c-b).length() ** 2)
        # if 2* area_tri < EPSILON * area_square: return True
        return b.minus(self).cross(c.minus(self)).is_zero()

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

    def is_strictly_within(self, p, r):
        """
        Return true if q is between p and r (exclusive).
        >>> Vector(1,0,0).is_within(Vector(1,0,0), Vector(2,0,0))
        True
        >>> Vector(1,0,0).is_strictly_within(Vector(1,0,0), Vector(2,0,0))
        False
        """
        return (p.x < self.x < r.x or r.x < self.x < p.x) and \
               (p.y < self.y < r.y or r.y < self.y < p.y) and \
               (p.z < self.z < r.z or r.z < self.z < p.z)


class Plane():
    def __init__(self, normal, distance):
        self.normal = Vector(normal).unit()
        self.distance = distance  # distance from (0, 0, 0)

    def __repr__(self):
        return 'Plane(normal={}, distance={})'.format(
                self.normal, self.distance
                )

    def clone(self):
        return Plane(self.normal.clone(), self.distance)

    @classmethod
    def from_points(cls, points):  # FIXME check flat face
        """
        Get a Plane from a list of points, after checking for collinearity.
        If collinear, return a None.
        >>> Plane.from_points(((0,0,0,),(1,0,0),(2,0,0),(3,0,1)))
        Plane(normal=Vector(0.000, -1.000, 0.000), distance=0)
        >>> Plane.from_points(((0,0,0,),(1,0,0),(2,0,0),(3,0,0)))
        """
        len_points = len(points)
        for i, a in enumerate(points):
            a = Vector(a)
            b = Vector(points[(i + 1) % len_points])
            c = Vector(points[(i + 2) % len_points])
            normal = b.minus(a).cross(c.minus(a))
            if not normal.is_zero():
                return Plane(normal.unit(), normal.dot(a))
        return None

    def flip(self):
        """
        Flip self normal.
        >>> p = Plane(normal=(1,0,0), distance=5); p.flip() ; p
        Plane(normal=Vector(-1.000, -0.000, -0.000), distance=-5)
        """
        self.normal = -self.normal
        self.distance = -self.distance


geometry = [None, None, None, None, None, None, ]


class Geom():
    def __init__(self, verts=None, polygons=None):
        if verts is None:
            verts = ()
        if polygons is None:
            polygons = ((), )
        try:
            # [1.,2.,3., 2.,3.,4., ...]
            self.verts = array.array('f', verts)
            # [[0,1,2,4], [1,2,3], ...]
            self.polygons = [list(p[:]) for p in polygons]
        except TypeError:
            raise Exception('Bad Geom() at init')

    def __repr__(self):
        """
        >>> Geom((1.,2.,3., 1.,2.,3.,), ((1,2,3),(1,2,3,4),(1,2,3,4,5),), )
        Geom(
            (1.000,2.000,3.000,  1.000,2.000,3.000),
            [[1, 2, 3], [1, 2, 3, 4], [1, 2, 3, 4, 5]],
            )
        """
        strverts = ('{:.3f}'.format(v) for v in self.verts)
        strverts = list(zip(*[iter(strverts)] * 3))  # Join: ((1.,2.,3.,), ...)
        strverts = ',  '.join((','.join(v) for v in strverts))
        return 'Geom(\n    ({}),\n    {},\n    )'.format(
                strverts, self.polygons,
                )

    def clone(self):
        return Geom(
                verts=self.verts[:],
                polygons=[p[:] for p in self.polygons],
                )

    def flip(self):
        """
        Flip self polygons normals.
        >>> g = Geom((), ((1,2,3),(1,2,3,4),(1,2,3,4,5),), ); g.flip(); g
        Geom(
            (),
            [[3, 2, 1], [4, 3, 2, 1], [5, 4, 3, 2, 1]],
            )
        """
        for p in self.polygons:
            p.reverse()

    def get_polygon(self, ipolygon):
        """
        Get ipolygon connectivity.
        >>> g = Geom((), ((1,2,3),(1,2,3,4),(1,2,3,4,5),), ); g.get_polygon(1)
        [1, 2, 3, 4]
        """
        return self.polygons[ipolygon]

    def get_polygon_verts(self, ipolygon):
        """
        Get ipolygon verts.
        >>> g = Geom((-1,-1,0, 1,-1,0, 0,1,0, 0,0,1), \
                     ((2,1,0), (0,1,3), (1,2,3), (2,0,3)) )  # Good tet
        >>> g.get_polygon_verts(0)
        [Vector(0.000, 1.000, 0.000), Vector(1.000, -1.000, 0.000), Vector(-1.000, -1.000, 0.000)]
        """
        return [self.get_vert(ivert) for ivert in self.get_polygon(ipolygon)]

    def update_polygon(self, ipolygon, polygon):
        """
        Update ipolygon connectivity.
        >>> g = Geom((-1,-1,0, 1,-1,0, 0,1,0, 0,0,1), \
                     ((2,1,0), (0,1,3), (1,2,3), (2,0,3)) )  # Good tet
        >>> g.update_polygon(0, (0,0,0)); g.get_polygon(0)
        0
        [0, 0, 0]
        """
        self.polygons[ipolygon] = list(polygon)
        return ipolygon

    def append_polygon(self, polygon):
        """
        Append a polygon.
        >>> g = Geom((-1,-1,0, 1,-1,0, 0,1,0, 0,0,1), \
                     ((2,1,0), (0,1,3), (1,2,3), (2,0,3)) )  # Good tet
        >>> g.append_polygon((0,0,0)); g.get_polygon(4)
        4
        [0, 0, 0]
        """
        self.polygons.append(list(polygon))
        ipolygon = self.get_npolygons() - 1
        return ipolygon

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

    def get_plane_of_polygon(self, ipolygon):
        """
        Get plane containing ipolygon.
        >>> g = Geom((0,0,1, 1,0,1, 2,0,1, 0,1,1,), ((0,1,2,3), ))
        >>> g.get_plane_of_polygon(0)
        Plane(normal=Vector(0.000, -0.000, 1.000), distance=1.0)
        """
        return Plane.from_points(self.get_polygon_verts(ipolygon))

    def split_polygon(self, ipolygon, plane, coplanar_front, coplanar_back, front, back):  # FIXME
        """
        Split ipolygon by a plane. Put the fragments in the inline lists.
        """
        COPLANAR = 0 # vertex of polygon within EPSILON distance from plane
        FRONT = 1    # vertex of polygon in front of the plane
        BACK = 2     # vertex of polygon at the back of the plane
        SPANNING = 3 # spanning polygon
        
        polygon = self.get_polygon(ipolygon)
        polygon_verts = self.get_polygon_verts(ipolygon)
        polygon_type = 0
        ivert_types = []       
        polygon_nverts = len(polygon)

        # Classify iverts types
        for i in range(polygon_nverts):
            distance = plane.normal.dot(polygon_verts[i]) - plane.distance
            ivert_type = -1
            if distance < -EPSILON: 
                ivert_type = BACK
            elif distance > EPSILON: 
                ivert_type = FRONT
            else: 
                ivert_type = COPLANAR
            polygon_type |= ivert_type
            ivert_types.append(ivert_type)
    
        # Classify the polygon and put it in the correct list
        if polygon_type == COPLANAR:
            polygon_plane = self.get_plane_of_polygon(ipolygon)          
            if plane.normal.dot(polygon_plane.normal) > 0:
                coplanar_front.append(polygon)
            else:
                coplanar_back.append(polygon)
        elif polygon_type == FRONT:
            front.append(polygon)
        elif polygon_type == BACK:
            back.append(polygon)
        elif polygon_type == SPANNING:  # FIXME FIXME arrivato fino a qui!
            f = []
            b = []
            for i in range(polygon_nverts):
                j = (i+1) % polygon_nverts
                ti = ivert_types[i]
                tj = ivert_types[j]
                vi = polygon_verts[i]
                vj = polygon_verts[j]
                if ti != BACK: 
                    f.append(vi)
                if ti != FRONT:
                    if ti != BACK: 
                        b.append(vi.clone())
                    else:
                        b.append(vi)
                if (ti | tj) == SPANNING:
                    # interpolation weight at the intersection point
                    t = (plane.distance - plane.normal.dot(vi)) / plane.normal.dot(vj.minus(vi))
                    # intersection point on the plane
                    v = vi.interpolate(vj, t)
                    f.append(v)
                    b.append(v.clone())
            if len(f) >= 3: 
                front.append(Polygon(f, polygon.shared))
            if len(b) >= 3: 
                back.append(Polygon(b, polygon.shared))




    def split_iface(igeom, iface, plane):
        """
        Split iface from igeom by spl_iface of spl_igeom if needed.
        Append new verts and new faces to geometry igeom.
        Return ifaces in the appropriate lists.
        >>> geometry[0] = Geom([-1,0,0, 1,0,0, 0,1,0, -1,1,0, 1,1,0],[0,1,2, 0,2,3, 1,4,2]) # z = 0, axis +k
        >>> split_iface(igeom=0, iface=0, plane=Plane((1,0,0),.5))  # other faces to cut
        ([], [], [0], [3, 4], {2: 5})
        >>> geometry[0]
        Geom(
             array('f', [-1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 0.0, 1.0, 1.0, 0.0, 0.5, 0.0, 0.0, 0.5, 0.5, 0.0]),
             array('i', [5, 1, 6, 0, 2, 3, 2, 6, 4, 0, 5, 6, 0, 6, 2, 6, 1, 4]),
        )
        >>> geometry[0] = Geom([-1,0,0, 1,0,0, 0,1,0, -1,1,0, 1,1,0],[0,1,2, 0,2,3, 1,4,2]) # axis +z
        >>> split_iface(igeom=0, iface=0, plane=Plane((1,0,0),-.5))  # other faces to cut
        ([], [], [0, 3], [4], {1: 5})
        >>> geometry[0]
        Geom(
             array('f', [-1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 0.0, 1.0, 1.0, 0.0, -0.5, 0.0, 0.0, -0.5, 0.5, 0.0]),
             array('i', [5, 1, 2, 0, 6, 3, 1, 4, 2, 5, 2, 6, 0, 5, 6, 6, 2, 3]),
        )
        >>> geometry[0] = Geom([-1,0,0, 1,0,0, 0,1,0, -1,1,0, 1,1,0],[0,1,2, 0,2,3, 1,4,2]) # axis +z
        >>> split_iface(igeom=0, iface=0, plane=Plane((0,-1,0),-.5))  # other faces to cut
        ([], [], [0, 3], [4], {2: 5, 1: 6})
        >>> geometry[0]
        Geom(
             array('f', [-1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 0.0, 1.0, 1.0, 0.0, 0.5, 0.5, 0.0, -0.5, 0.5, 0.0]),
             array('i', [0, 1, 5, 0, 6, 3, 2, 5, 4, 0, 5, 6, 5, 2, 6, 5, 1, 4, 6, 2, 3]),
        )
        >>> geometry[0] = Geom([-1,0,0, 1,0,0, 0,1,0],[0,1,2])  # axis +k
        >>> split_iface(igeom=0, iface=0, plane=Plane((0,0,-1),0))  # copl back
        ([], [0], [], [], {})
        >>> split_iface(igeom=0, iface=0, plane=Plane((0,0, 1),0))  # copl front
        ([0], [], [], [], {})
        """
        # Vertices and faces types, collections of ifaces
        COPLANAR = 0  # vertex or face is within EPSILON2 distance from plane
        FRONT = 1     # vertex or face is in front of the plane
        BACK = 2      # vertex or face is at the back of the plane
        SPANNING = 3  # edge is intersected
        coplanar_front, coplanar_back, front, back = [], [], [], []
    
        # The edges to be split,
        # eg. {(2,3): 1} with {(ivert0,ivert1): cut_ivert, ...}
        # The opposite of the current edge is sent for easy search
        spl_edges = {}
    
        # Classify each point as well as the entire polygon
        # into one of the above classes.
        faceType = 0
        vertexTypes = []
        spl_normal, spl_distance = plane.n, plane.w
        face = get_face(igeom, iface)
        for ivert in face:
            # Calc the distance between the vertex and the splitting plane
            # then classify the vertex, and update classification of the face
            vert = get_vert(igeom, ivert)
            distance = spl_normal.dot(vert) - spl_distance
            if distance < -EPSILON2:
                vertexType = BACK
            elif distance > EPSILON2:
                vertexType = FRONT
            else:
                vertexType = COPLANAR
            faceType |= vertexType
            vertexTypes.append(vertexType)
    
        # Put the face in the correct list, splitting it when necessary.
        if faceType == COPLANAR:
            iface_plane = get_plane_from_iface(igeom, iface)
            iface_normal = iface_plane.n
            if spl_normal.dot(iface_normal) > 0:
                coplanar_front.append(iface)
            else:
                coplanar_back.append(iface)
        elif faceType == FRONT:
            front.append(iface)
        elif faceType == BACK:
            back.append(iface)
        elif faceType == SPANNING:  # the face is spanning, so cut it
            front_iverts, back_iverts = [], []  # front vertices, back vertices
            for i in range(3):
                j = (i+1) % 3
                ti = vertexTypes[i]  # can be BACK, COPLANAR, FRONT
                tj = vertexTypes[j]
                vi = face[i]
                vj = face[j]
                # the edge is on one side or spanning?
                if ti != BACK:
                    front_iverts.append(vi)
                if ti != FRONT:
                    if ti != BACK:
                        back_iverts.append(vi)
                    else:
                        back_iverts.append(vi)
                if (ti | tj) == SPANNING:
                    # interpolation weight at the intersection point
                    vi_vert = get_vert(igeom, vi)
                    vj_vert = get_vert(igeom, vj)
                    t = (spl_distance - spl_normal.dot(vi_vert)) \
                        / spl_normal.dot(vj_vert.minus(vi_vert))
                    # intersection point on the plane
                    cut_vert = vi_vert.lerp(vj_vert, t)
                    cut_ivert = append_vert(igeom, cut_vert)
                    front_iverts.append(cut_ivert)
                    back_iverts.append(cut_ivert)
                    # update spl_edges, the opposite for later search!
                    spl_edges[(vj, vi)] = cut_ivert
    
            # Add front cut faces
            updated = False
    
            if len(front_iverts) > 2:
                if updated:
                    new_iface = append_face(
                        igeom,
                        front_iverts[0:3],
                        )
                else:
                    updated = True
                    new_iface = update_face(
                        igeom,
                        iface,
                        front_iverts[0:3],
                        )
                front.append(new_iface)
    
                if len(front_iverts) == 4:
                    if updated:
                        new_iface = append_face(
                            igeom,
                            (front_iverts[0], front_iverts[2], front_iverts[3]),
                            )
                    else:
                        updated = True
                        new_iface = update_face(
                            igeom,
                            iface,
                            (front_iverts[0], front_iverts[2], front_iverts[3]),
                            )
                    front.append(new_iface)
            else:
                raise Exception('Problem with front spanning:', front_iverts)
    
            # Add back cut faces
            if len(back_iverts) > 2:
                if updated:
                    new_iface = append_face(
                        igeom,
                        back_iverts[0:3],
                        )
                else:
                    updated = True
                    new_iface = update_face(
                        igeom,
                        iface,
                        back_iverts[0:3],
                        )
                back.append(new_iface)
                if len(back_iverts) == 4:
                    if updated:
                        new_iface = append_face(
                            igeom,
                            (back_iverts[0], back_iverts[2], back_iverts[3]),
                            )
                    else:
                        updated = True
                        new_iface = update_face(
                            igeom,
                            iface,
                            (back_iverts[0], back_iverts[2], back_iverts[3]),
                            )
                    back.append(new_iface)
            else:
                raise Exception('Problem with back spanning:', back_iverts)
    
        # Manage spl_edges FIXME FIXME test
        spl_ifaces = {}  # {iface: iface0, iface1}
        if spl_edges:  # FIXME
            halfedges = get_halfedges(igeom, get_ifaces(igeom))  # FIXME not all ifaces!
            for spl_edge, spl_ivert in spl_edges.items():
                # Get the bordering face that is split
                spl_iface = halfedges.get(spl_edge, None)
                if spl_iface is None:  # FIXME on a border
                    continue
                # Get its iverts, and find the far ivert
                spl_face = get_face(igeom, spl_iface)
                spl_face.remove(spl_edge[0])  # Remove unwanted
                spl_face.remove(spl_edge[1])  # Remove unwanted
                spl_face_ivert2 = spl_face[0]  # Get the left one
                # Build the two new faces
                new_iface0 = update_face(
                        igeom,
                        spl_iface,
                        (spl_edge[0], spl_ivert, spl_face_ivert2),
                    )
                new_iface1 = append_face(
                        igeom,
                        (spl_ivert, spl_edge[1], spl_face_ivert2),
                    )
                # Tell the caller about spl_iface1 (Only the new one!)
                spl_ifaces[spl_iface] = new_iface1
                
        # Return
        return coplanar_front, coplanar_back, front, back, spl_ifaces



















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


    def merge_duplicated_verts(self):
        """
        Remove dup verts, and relink all polygons. No mod of polygons.
        >>> g = Geom((-1,-1,0, 1,-1,0, 0,1,0, 0,0,1, 0,1,0, 0,1,0, 1,-1,0, 1,-1,0,), \
                     ((2,6,0), (0,1,3), (7,4,3), (5,0,3)) )  # Dup verts, 8>4
        >>> g.merge_duplicated_verts(); g
        Geom(
            (-1.000,-1.000,0.000,  1.000,-1.000,0.000,  0.000,1.000,0.000,  0.000,0.000,1.000),
            [[2, 1, 0], [0, 1, 3], [1, 2, 3], [2, 0, 3]],
            )
        """
        # Find unique verts and build ivert_to_ivert dict
        unique_pyverts = []
        ivert_to_ivert = {}
        nverts = self.get_nverts()
        for ivert in range(nverts):
            vert = self.get_vert(ivert)
            seen = False
            for i, selected_pyvert in enumerate(unique_pyverts):
                if (selected_pyvert - vert).is_zero():
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

    # edges

    def get_halfedges(self, ipolygons=None):
        """
        Get halfedges dict of ipolygons subset.
        halfedges are: {(1,2):7]} with {(ivert0, ivert1): ipolygon on the left}
        according to iface0 normal up
        >>> g = Geom((-1,-1,0, 1,-1,0, 0,1,0, 0,0,1), \
                     ((2,1,0), (0,1,3), (1,2,3), (2,0,3)) )  # Good tet
        >>> g.get_halfedges(ipolygons=None)
        {(2, 1): 0, (1, 0): 0, (0, 2): 0, (0, 1): 1, (1, 3): 1, (3, 0): 1, (1, 2): 2, (2, 3): 2, (3, 1): 2, (2, 0): 3, (0, 3): 3, (3, 2): 3}
        >>> g.get_halfedges(ipolygons=(1,2,3))
        {(0, 1): 1, (1, 3): 1, (3, 0): 1, (1, 2): 2, (2, 3): 2, (3, 1): 2, (2, 0): 3, (0, 3): 3, (3, 2): 3}
        >>> g = Geom((-1,-1,0, 1,-1,0, 0,1,0, 0,0,1), \
                     ((2,1,0), (0,1,3), (1,2,3), (3,0,2)) )  # Unorient tet
        >>> g.get_halfedges()
        Traceback (most recent call last):
        ...
        Exception: ('Invalid GEOM, non-manifold or unorientable at ipolygon:', 3)
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
                    raise Exception('Invalid GEOM, non-manifold or unorientable at ipolygon:', ipolygon)
                halfedges[halfedge] = ipolygon
        return halfedges

    def get_border_halfedges(self, ipolygons=None):
        """
        Get border halfedges dict
        Eg: {(1,2):7]} with {(ivert0, ivert1): iface on the left}
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

    def get_tris_of_polygon(self, ipolygon):
        """
        Triangulate ipolygon with no zero-area tris
        >>> g = Geom((0,0,0, 1,0,0, 2,0,0, 3,0,0, 1,1,0, 0,1,0), \
                     ((0,1,2,3,4,5), ))  # Polyhedra, 6 edges, 3 collinear
        >>> g.get_tris_of_polygon(ipolygon=0)
        [(4, 5, 0), (4, 0, 1), (4, 1, 2), (4, 2, 3)]
        >>> g = Geom((0,0,0, 1,0,0, 2,0,0, 3,0,0), \
                     ((0,1,2,3,), ))    # Zero area polyhedra
        >>> g.get_tris_of_polygon(ipolygon=0)
        Traceback (most recent call last):
        ...
        Exception: ('Zero area or edge, ipolygon:', 0)
        >>> g = Geom((0,0,0, 1,0,0, 1,0,0, 3,1,0), \
                     ((0,1,2,3,), ))    # Zero lenght edge polyhedra
        >>> g.get_tris_of_polygon(ipolygon=0)
        Traceback (most recent call last):
        ...
        Exception: ('Zero area or edge, ipolygon:', 0)
        """
        polygon = self.get_polygon(ipolygon)
        polygon_nverts = len(polygon)
        choices = {}
        for start in range(polygon_nverts):
            tot_perimeter = 0.
            zero_area = False
            tris = []
            for i in range(polygon_nverts-2):
                # Build vertices
                ia = polygon[start % polygon_nverts]
                ib = polygon[(start + i + 1) % polygon_nverts]
                ic = polygon[(start + i + 2) % polygon_nverts]
                a = self.get_vert(ia)
                b = self.get_vert(ib)
                c = self.get_vert(ic)
                # Check not zero area
                cross = b.minus(a).cross(c.minus(a))
                if cross.is_zero():
                    zero_area = True
                    break
                tot_perimeter += a.length() + b.length() + c.length()
                tris.append((ia, ib, ic))
            if not zero_area:
                choices[tot_perimeter] = tris
        if choices:
            return choices[min(choices)]
        else:
            raise Exception('Zero area or edge, ipolygon:', ipolygon)

    # STL

    def to_STL(self, filename):
        """
        Write self to STL file
        >>> g = Geom((-1.0, -1.0, -1.0,  -1.0, -1.0, 1.0,  -1.0, 1.0,  1.0, \
                      -1.0,  1.0, -1.0,   1.0,  1.0, 1.0,   1.0, 1.0, -1.0,\
                       1.0, -1.0, -1.0,   1.0, -1.0, 1.0), \
                    ((0,1,2,3), (7,6,5,4), (1,7,4,2), \
                     (0,3,5,6), (1,0,6,7), (2,4,5,3)) )  # A good cube
        >>> g.to_STL('../test/doctest.stl')
        to_STL: ../test/doctest.stl
        """
        with open(filename, 'w') as f:
            f.write('solid name\n')
            for ipolygon in range(self.get_npolygons()):
                tris = self.get_tris_of_polygon(ipolygon)
                for tri in tris:
                    f.write('facet normal 0 0 0\n')
                    f.write(' outer loop\n')
                    for ivert in tri:
                        f.write('  vertex {v[0]:.9f} {v[1]:.9f} {v[2]:.9f}\n'.format(
                                v=self.get_vert(ivert)))
                    f.write(' endloop\n')
                    f.write('endfacet\n')
            f.write('endsolid name\n')
        print('to_STL:', filename)

    @classmethod
    def from_STL(cls, filename):
        """
        Get new Geom from STL file
        >>> Geom.from_STL('../test/doctest.stl')
        Geom(
            (-1.000,1.000,-1.000,  -1.000,-1.000,-1.000,  -1.000,-1.000,1.000,  -1.000,1.000,1.000,  1.000,1.000,1.000,  1.000,-1.000,1.000,  1.000,-1.000,-1.000,  1.000,1.000,-1.000),
            [[0, 1, 2], [0, 2, 3], [4, 5, 6], [4, 6, 7], [3, 2, 5], [3, 5, 4], [6, 1, 0], [6, 0, 7], [5, 2, 1], [5, 1, 6], [0, 3, 4], [0, 4, 7]],
            )
        """
        # Get STL mesh
        from stl import mesh
        mesh = mesh.Mesh.from_file(filename)
        verts, polygons, py_verts = [], [], []
        for iface, p in enumerate(mesh.points):
            # p is [-1.0, -1.0, -1.0, -1.0, -1.0, 1.0, -1.0, 1.0, 1.0]
            verts.extend(p)
            polygons.append((3*iface, 3*iface+1, 3*iface+2))
            py_verts.append((p[0], p[1], p[2]))
            py_verts.append((p[3], p[4], p[5]))
            py_verts.append((p[6], p[7], p[8]))
        g = Geom(verts, polygons)
        g.merge_duplicated_verts()
        g.check_geom_sanity()
        return g

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
        for ipolygon in range(self.get_npolygons()):
            self.get_tris_of_polygon(ipolygon)

    def check_is_solid(self):
        """
        Check surface:
        - 2-manifold and closed, each edge should join two faces, no more no less
        - orientable, adjoining faces should have normals in the same directions
        >>> g = Geom((-1,-1,0, 1,-1,0, 0,1,0, 0,0,1), \
                     ((0,1,3), (1,2,3), (2,0,3)) )  # Open tet
        >>> g.check_is_solid()
        Traceback (most recent call last):
        ...
        Exception: ('Invalid GEOM, non closed at polygons:', [2, 1, 0])
        """
        border_halfedges = self.get_border_halfedges()
        if border_halfedges:
            raise Exception("Invalid GEOM, non closed at polygons:",
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
#
#
    def check_geom_sanity(self):
        """
        Check geometry sanity
    
        If the mesh is correct and encloses a volume, this can be checked with
        prior tests: checking orientability, non-borders, non-self-intersecting.
        After that we can calculate its topological features and check if
        Euler's formula  C+V=A+2(S-H) is satisfied.
        If the mesh is not correct, many geometric algorithms will fail.
        The only solution in this case is the user repairing the mesh.
        >>> g = Geom((-1.0, -1.0, -1.0,  -1.0, -1.0, 1.0,  -1.0, 1.0,  1.0, \
                      -1.0,  1.0, -1.0,   1.0,  1.0, 1.0,   1.0, 1.0, -1.0,\
                       1.0, -1.0, -1.0,   1.0, -1.0, 1.0), \
                    ((0,1,2,3), (7,6,5,4), (1,7,4,2), \
                     (0,3,5,6), (1,0,6,7), (2,4,5,3)) )  # A good cube
        >>> g.check_geom_sanity()
        """
        self.check_loose_verts()
        self.check_degenerate_geometry()
        self.check_is_solid()
#        self.check_euler()  
        # Check correct normals for a solid in fluid FIXME working here
        # Check self intersection  # FIXME not working
    























if __name__ == "__main__":
    import doctest
    doctest.testmod()
