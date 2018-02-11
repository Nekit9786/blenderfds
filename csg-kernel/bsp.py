#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 10 17:42:57 2018

@author: egissi
"""
import textwrap
from geometry import *


class BSP():
    """
    BSP tree of a geometry
    """
    def __init__(self, igeom, ifaces=None):
        self.igeom = igeom     # Referred igeom
        self.spl_iface = None  # Reference to splitting iface
        self.ifaces = []       # Coplanar ifaces
        self.front_bsp = None  # link to child front BSP
        self.back_bsp = None   # link to child back BSP
        if ifaces:
            self.build(ifaces)

    def __repr__(self):
        return 'BSP tree of{}'.format(self._repr_tree())

    def _repr_tree(self):
        # Get children trees
        front_tree = "None"
        if self.front_bsp:
            front_tree = self.front_bsp._repr_tree()
        back_tree = "None"
        if self.back_bsp:
            back_tree = self.back_bsp._repr_tree()
        # Join texts
        text = 'igeom: {}, spl_iface: {}, ifaces: {}\n└─front_bsp: {}\n└─back_bsp:  {}'.format(
            self.igeom,
            self.spl_iface,
            self.ifaces,
            front_tree,
            back_tree,
        )
        return textwrap.indent(text, '  ')

    def clone(self):
        bsp = BSP(igeom=self.igeom)
        bsp.spl_iface = self.spl_iface
        bsp.ifaces = self.ifaces[:]
        if self.front_bsp:
            bsp.front_bsp = self.front_bsp.clone()
        if self.back_bsp:
            bsp.back_bsp = self.back_bsp.clone()
        return bsp

    def build(self, ifaces):
        """
        Build self from ifaces
        >>> geometry[0] = Geom([-1.0, -1.0, -1.0, -1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0, 1.0, -1.0, 1.0, 1.0, 1.0, 1.0, 1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0, -1.0], [0, 1, 2, 2, 3, 0, 3, 2, 4, 4, 5, 3, 5, 4, 6, 6, 7, 5, 1, 0, 7, 7, 6, 1, 7, 0, 3, 3, 5, 7, 4, 2, 1, 1, 6, 4])  # A cube
        >>> print(BSP(igeom=0, ifaces=get_ifaces(0)))
        BSP tree of  igeom: 0, spl_iface: 0, ifaces: [0]
          └─front_bsp:   igeom: 0, spl_iface: 1, ifaces: [1]
            └─front_bsp: None
            └─back_bsp:  None
          └─back_bsp:    igeom: 0, spl_iface: 2, ifaces: [2]
            └─front_bsp:   igeom: 0, spl_iface: 3, ifaces: [3]
              └─front_bsp: None
              └─back_bsp:  None
            └─back_bsp:    igeom: 0, spl_iface: 4, ifaces: [4]
              └─front_bsp:   igeom: 0, spl_iface: 5, ifaces: [5]
                └─front_bsp: None
                └─back_bsp:  None
              └─back_bsp:    igeom: 0, spl_iface: 6, ifaces: [6]
                └─front_bsp:   igeom: 0, spl_iface: 7, ifaces: [7]
                  └─front_bsp: None
                  └─back_bsp:  None
                └─back_bsp:    igeom: 0, spl_iface: 8, ifaces: [8]
                  └─front_bsp:   igeom: 0, spl_iface: 9, ifaces: [9]
                    └─front_bsp: None
                    └─back_bsp:  None
                  └─back_bsp:    igeom: 0, spl_iface: 10, ifaces: [10]
                    └─front_bsp:   igeom: 0, spl_iface: 11, ifaces: [11]
                      └─front_bsp: None
                      └─back_bsp:  None
                    └─back_bsp:  None
        """
        # Protect
        if not ifaces:
            return
        # Use first iface as splitting iface  # FIXME check!
        i = 0
        if not self.spl_iface:
            self.spl_iface = ifaces[0]
            self.ifaces.append(ifaces[0])
            i = 1
        # Select ifaces for front and back, split them if needed.
        # front and back are lists of ifaces
        front = []
        back = []
        for iface in ifaces[i:]:
            # coplanar front and back polygons go into self.polygons  # FIXME new_cut_iverts?
            new_coplanar_front, new_coplanar_back, new_front, new_back, new_cut_iverts = split_iface(
                    igeom=self.igeom,
                    iface=iface,
                    spl_igeom=self.igeom,
                    spl_iface=self.spl_iface,
                    )
            front.extend(new_front)
            back.extend(new_back)
            front.extend(new_coplanar_front)
            back.extend(new_coplanar_back)
        # Recursively build the BSP tree and return it
        if len(front) > 0:
            if not self.front_bsp:
                self.front_bsp = BSP(self.igeom)
            self.front_bsp.build(ifaces=front)
        if len(back) > 0:
            if not self.back_bsp:
                self.back_bsp = BSP(self.igeom)
            self.back_bsp.build(ifaces=back)

    def invert(self):
        """
        Convert self solid space to empty space and empty space to solid space.
        >>> geometry[0] = Geom([-1.0, -1.0, -1.0, -1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0, 1.0, -1.0, 1.0, 1.0, 1.0, 1.0, 1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0, -1.0], [0, 1, 2, 2, 3, 0, 3, 2, 4, 4, 5, 3, 5, 4, 6, 6, 7, 5, 1, 0, 7, 7, 6, 1, 7, 0, 3, 3, 5, 7, 4, 2, 1, 1, 6, 4])  # A cube
        >>> BSP(igeom=0, ifaces=get_ifaces(0)).invert()
        >>> print(get_iface_plane(0, 5))
        (Vector(-1.000, 0.000, 0.000), -1.0)
        """
        igeom = self.igeom
        # Flip all face normals
        if self.spl_iface not in self.ifaces:
            # FIXME otherwise risk of double flip
            flip_iface_normal(igeom, self.spl_iface)
        for iface in self.ifaces:
            flip_iface_normal(igeom, iface)
        if self.front_bsp:
            self.front_bsp.invert()
        if self.back_bsp:
            self.back_bsp.invert()
        # Swap front and back bsp
        tmp = self.front_bsp
        self.front_bsp = self.back_bsp
        self.back_bsp = tmp


def clip_ifaces(igeom, ifaces, clipping_bsp):  # FIXME test
    """
    Recursively remove all ifaces that are inside the clipping_bsp tree.
    """
    if clipping_bsp.spl_iface is None:  # FIXME why? It always have it
        return ifaces[:]

    front = []
    back = []
    for iface in ifaces:
        new_coplanar_front, new_coplanar_back, new_front, new_back, new_cut_iverts = split_iface(
                igeom=igeom,
                iface=iface,
                spl_igeom=clipping_bsp.igeom,  # FIXME move out of for cycle
                spl_iface=clipping_bsp.spl_iface,    # FIXME move out of for cycle
                )
        front.extend(new_front)
        back.extend(new_back)
        front.extend(new_coplanar_front)
        back.extend(new_coplanar_back)

    if clipping_bsp.front_bsp:
        # recurse on branches, conserve those polygons
        front = clip_ifaces(igeom, front, clipping_bsp.front_bsp)

    if clipping_bsp.back_bsp:
        # recurse on branches, conserve those polygons
        back = clip_ifaces(igeom, back, clipping_bsp.back_bsp)
    else:
        back = []  # but remove polygons that are back of the leaves

    # send all non removed polygons, for recursion
    front.extend(back)
    return front


def clip_to(bsp, clipping_bsp):  # FIX name, it returns None   # FIXME test
    """
    Remove all polygons in bsp tree that are inside the clipping_bsp tree.
    Just send all bsp ifaces to clip_ifaces().
    """

    bsp.ifaces = clip_ifaces(
            igeom=bsp.igeom,
            ifaces=bsp.ifaces,
            clipping_bsp=clipping_bsp,
            )

    if bsp.front_bsp:
        clip_to(
                bsp=bsp.front_bsp,
                clipping_bsp=clipping_bsp,
                )
    if bsp.back_bsp:
        clip_to(
                bsp=bsp.back_bsp,
                clipping_bsp=clipping_bsp,
                )


def get_all_ifaces_from_bsp(bsp):
    """
    Get recursively all ifaces from bsp and its children
    >>> geometry[0] = Geom([-1.0, -1.0, -1.0, -1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0, 1.0, -1.0, 1.0, 1.0, 1.0, 1.0, 1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0, -1.0], [0, 1, 2, 2, 3, 0, 3, 2, 4, 4, 5, 3, 5, 4, 6, 6, 7, 5, 1, 0, 7, 7, 6, 1, 7, 0, 3, 3, 5, 7, 4, 2, 1, 1, 6, 4])  # A cube
    >>> get_all_ifaces_from_bsp(BSP(igeom=0, ifaces=get_ifaces(0)))
    [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
    """
    ifaces = bsp.ifaces[:]
    if bsp.front_bsp:
        ifaces.extend(get_all_ifaces_from_bsp(bsp.front_bsp))
    if bsp.back_bsp:
        ifaces.extend(get_all_ifaces_from_bsp(bsp.back_bsp))
    ifaces.sort()
    return ifaces


# Boolean operations


def geom_union(igeom0, igeom1, name='union'):  # FIXME
    """
    Union of igeom0 and igeom1
#    >>> geometry[0] = Geom([-1,-1,0,   1,-1,0,   0,1,0,   0,0,1],   [2,1,0, 0,1,3, 1,2,3, 2,0,3])  # Good tet
#    >>> geometry[1] = Geom([-.5,-1,0,  1.5,-1,0, 0.5,1,0, 0.5,0,1], [2,1,0, 0,1,3, 1,2,3, 2,0,3])  # Good tet, move x .5
#    >>> geom_union(0, 1)
#    Geom(
#         array('f', [-1.0, -1.0, 0.0, 1.0, -1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.25, 0.5, 0.0, 0.25, -0.25, 0.75, -0.5, -1.0, 0.0, 1.5, -1.0, 0.0, 0.5, 1.0, 0.0, 0.5, 0.0, 1.0, 0.25, 0.5, 0.0, 0.25, -0.25, 0.75, 1.0, -1.0, 0.0, 0.25, 0.5, 0.0, 1.0, -1.0, 0.0, 0.25, -0.25, 0.75]),
#         array('i', [2, 1, 0, 0, 1, 3, 2, 0, 3, 4, 2, 3, 4, 3, 5, 7, 8, 9, 8, 10, 11, 8, 11, 9, 8, 7, 12, 13, 8, 12, 14, 7, 9, 15, 14, 9]),
#    )
    """
    a = BSP(igeom=igeom0, ifaces=get_ifaces(igeom0))
    b = BSP(igeom=igeom1, ifaces=get_ifaces(igeom1))

    # Clip
    clip_to(a, b)  # remove everything in a inside b
    clip_to(b, a)  # remove everything in b inside a
    b.invert()
    clip_to(b, a)  # remove everything in -b inside a
    b.invert()

    # Create new geometry
    geometry[igeom0] = get_new_geom_from_bsp(a)
    geometry[igeom1] = get_new_geom_from_bsp(b)

    # check_geom_sanity(igeom0)

    # Send to STL
    to_STL(igeom0, filename='{}_a_clipped.stl'.format(name))
    to_STL(igeom1, filename='{}_b_clipped.stl'.format(name))

    # Join
#    geometry[igeom0] = get_joined_geom(igeom0, igeom1)

    # Send to STL
#    to_STL(igeom0, filename='{}_ab.stl'.format(name))

    return geometry[igeom0]

#def geom_intersection(igeom0, igeom1, name='inters'):  # FIXME working on this
#    """
#    Intersection of igeom0 and igeom1
#    >>> geometry[0] = Geom([-1,-1,0,   1,-1,0,   0,1,0,   0,0,1],   [2,1,0, 0,1,3, 1,2,3, 2,0,3])  # Good tet
#    >>> geometry[1] = Geom([-.5,-1,0,  1.5,-1,0, 0.5,1,0, 0.5,0,1], [2,1,0, 0,1,3, 1,2,3, 2,0,3])  # Good tet, move x .5
#    >>> geom_intersection(0, 1)
#    """
#    a = BSP(igeom=igeom0, ifaces=get_ifaces(igeom0))
#    b = BSP(igeom=igeom1, ifaces=get_ifaces(igeom1))
#
#    # Send to STL
#    to_STL(igeom0, filename='{}_a.stl'.format(name))
#    to_STL(igeom1, filename='{}_b.stl'.format(name))
#
#    a.invert()
#    clip_to(b, a)  # remove everything in b inside -a
#    b.invert()
#    clip_to(a, b)  # remove everything in -a inside -b
#    clip_to(b, a)  # remove everything in -b inside -a
#
#    a.invert()
#    b.invert()
#
#    # Create new geometry
#    geometry[igeom0] = get_new_geom_from_bsp(a)
#    geometry[igeom1] = get_new_geom_from_bsp(b)
#
#    # Send to STL
#    to_STL(igeom0, filename='{}_a_clipped.stl'.format(name))
#    to_STL(igeom1, filename='{}_b_clipped.stl'.format(name))
#
#    # Join
#    geometry[igeom0] = get_joined_geom(igeom0, igeom1)
#
#    # Send to STL
#    to_STL(igeom0, filename='{}_ab.stl'.format(name))
#
#    return geometry[igeom0]

