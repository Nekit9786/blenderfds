# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 08:42:30 2018

@author: eng4
"""

import textwrap
from geometric_entities import Geom


class BSPNode(object):
    """
    A node in a BSP tree.
    >>> g = Geom((-1,-1,0, 1,-1,0, 0,1,0, 0,0,1), \
                 ((2,1,0), (0,1,3), (1,2,3), (2,0,3)) )  # Good tet
    >>> n = BSPNode(geom=g); n.build(); n
    BSP tree - geom hid: None, ipolygons: [0]
        │ plane: Plane(normal=Vector(0.000, -0.000, -1.000), distance=0.0)
        ├─front_node: None
        └─back_node: geom hid: None, ipolygons: [1]
          │ plane: Plane(normal=Vector(0.000, -0.707, 0.707), distance=2.0)
          ├─front_node: None
          └─back_node: geom hid: None, ipolygons: [2]
            │ plane: Plane(normal=Vector(0.816, 0.408, 0.408), distance=1.0)
            ├─front_node: None
            └─back_node: geom hid: None, ipolygons: [3]
              │ plane: Plane(normal=Vector(-0.816, 0.408, 0.408), distance=1.0)
              ├─front_node: None
              └─back_node: None
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
        text = '{3}│ plane: {0}\n{3}├─front_node: {1}\n{3}└─back_node: {2}'.format(
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

    def invert(self):
        """
        Swap solid space and empty space.
        """
        # Flip normals
        self.geom.flip()
        self.plane.flip()
        # Invert tree
        if self.front_node:
            self.front_node.invert()
        if self.back_node:
            self.back_node.invert()
        # Swap front and back nodes
        temp = self.front_node
        self.front_node = self.back_node
        self.back_node = temp

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
            ipolygons.extend(self.front_node.all_ipolygons())
        if self.back_node:
            ipolygons.extend(self.back_node.all_ipolygons())
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
            self.plane = self.geom.get_plane_of_polygon(ipolygons[0])
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

    def append(self, other_bsp):  # FIXME test
        """
        Append other bsp tree to this.
        """
        # Append the other geom to self.geom
        new_ipolygons = self.geom.append(other_bsp.geom)
        # Build
        self.build(new_ipolygons)

if __name__ == "__main__":
    import doctest
    doctest.testmod()
