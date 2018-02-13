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
    def __init__(self, igeom, ifaces=None, master_bsp=None):
        self.igeom = igeom     # Referred igeom
        self.plane = None      # Reference to splitting plane
        self.ifaces = []       # Coplanar ifaces
        self.front_bsp = None  # link to child front BSP
        self.back_bsp = None   # link to child back BSP
        if master_bsp:
            self.master_bsp = master_bsp
        else:
            self.master_bsp = self
        if ifaces:
            self.build(ifaces)
        geometry[igeom].master_bsp = self.master_bsp  # FIXME this is an hack

    def __repr__(self):
        return 'BSP tree {}'.format(self._repr_tree())

    def _repr_tree(self, back=True):
        # Get children trees
        front_tree = "None"
        if self.front_bsp:
            front_tree = self.front_bsp._repr_tree(back=False)
        back_tree = "None"
        if self.back_bsp:
            back_tree = self.back_bsp._repr_tree(back=True)
        # Join texts
        line = "│ "
        if back:
            line = "  "
        header = 'igeom: {0}, ifaces: {1}\n'.format(
            self.igeom,
            self.ifaces,
            )
        text = '{5}│ plane: {1}\n{5}├─front_bsp {3}\n{5}└─back_bsp {4}'.format(
            self.igeom,
            self.plane,
            self.ifaces,
            front_tree,
            back_tree,
            line,
        )
        return header + textwrap.indent(text, '  ')

    def clone(self):
        bsp = BSP(igeom=self.igeom, master_bsp=self.master_bsp)
        bsp.plane = self.plane
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
        >>> BSP(igeom=0, ifaces=get_ifaces(0))
        BSP tree igeom: 0, ifaces: [0]
            │ plane: Plane(n=Vector(-1.000, 0.000, 0.000), w=1.0)
            ├─front_bsp igeom: 0, ifaces: [1]
            │ │ plane: Plane(n=Vector(-1.000, 0.000, -0.000), w=1.0)
            │ ├─front_bsp None
            │ └─back_bsp None
            └─back_bsp igeom: 0, ifaces: [2]
              │ plane: Plane(n=Vector(0.000, 1.000, 0.000), w=1.0)
              ├─front_bsp igeom: 0, ifaces: [3]
              │ │ plane: Plane(n=Vector(0.000, 1.000, 0.000), w=1.0)
              │ ├─front_bsp None
              │ └─back_bsp None
              └─back_bsp igeom: 0, ifaces: [4]
                │ plane: Plane(n=Vector(1.000, 0.000, -0.000), w=1.0)
                ├─front_bsp igeom: 0, ifaces: [5]
                │ │ plane: Plane(n=Vector(1.000, 0.000, 0.000), w=1.0)
                │ ├─front_bsp None
                │ └─back_bsp None
                └─back_bsp igeom: 0, ifaces: [6]
                  │ plane: Plane(n=Vector(0.000, -1.000, 0.000), w=1.0)
                  ├─front_bsp igeom: 0, ifaces: [7]
                  │ │ plane: Plane(n=Vector(0.000, -1.000, 0.000), w=1.0)
                  │ ├─front_bsp None
                  │ └─back_bsp None
                  └─back_bsp igeom: 0, ifaces: [8]
                    │ plane: Plane(n=Vector(0.000, 0.000, -1.000), w=1.0)
                    ├─front_bsp igeom: 0, ifaces: [9]
                    │ │ plane: Plane(n=Vector(0.000, 0.000, -1.000), w=1.0)
                    │ ├─front_bsp None
                    │ └─back_bsp None
                    └─back_bsp igeom: 0, ifaces: [10]
                      │ plane: Plane(n=Vector(0.000, 0.000, 1.000), w=1.0)
                      ├─front_bsp igeom: 0, ifaces: [11]
                      │ │ plane: Plane(n=Vector(0.000, 0.000, 1.000), w=1.0)
                      │ ├─front_bsp None
                      │ └─back_bsp None
                      └─back_bsp None
        """
        # Protect
        if not ifaces:
            return
        # Use first iface as splitting iface  # FIXME check!
        i = 0
        if not self.plane:
            self.plane = get_plane_from_iface(self.igeom, ifaces[0])
            self.ifaces.append(ifaces[0])
            i = 1
        # Select ifaces for front and back, split them if needed.
        # front and back are lists of ifaces
        master_bsp = self.master_bsp
        front = []
        back = []
        for iface in ifaces[i:]:
            # coplanar front and back polygons go into self.polygons  # FIXME new_cut_iverts?
            new_coplanar_front, new_coplanar_back, new_front, new_back, spl_ifaces = split_iface(
                    igeom=self.igeom,
                    iface=iface,
                    plane=self.plane,
                    )
            front.extend(new_front)
            back.extend(new_back)
            front.extend(new_coplanar_front)
            back.extend(new_coplanar_back)
            
            # Update lists and bsp tree
            for spl_iface, new_iface in spl_ifaces.items():
                if spl_iface in ifaces:
                    ifaces.append(new_iface)
                elif spl_iface in front:
                    front.append(new_iface)
                elif spl_iface in back:
                    back.append(new_iface)                
                else:
                    ret = update_iface_in_bsp(master_bsp, spl_iface, new_iface)
                    if not ret:
                        print("B: NOT FOUND iface:", spl_iface)
        
        # Recursively build the BSP tree and return it
        if len(front) > 0:
            if not self.front_bsp:
                self.front_bsp = BSP(self.igeom, master_bsp=self.master_bsp)
            self.front_bsp.build(ifaces=front)
        if len(back) > 0:
            if not self.back_bsp:
                self.back_bsp = BSP(self.igeom, master_bsp=self.master_bsp)
            self.back_bsp.build(ifaces=back)

    def invert(self):
        """
        Convert self solid space to empty space and empty space to solid space.
        >>> geometry[0] = Geom([-1.0, -1.0, -1.0, -1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0, 1.0, -1.0, 1.0, 1.0, 1.0, 1.0, 1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0, -1.0], [0, 1, 2, 2, 3, 0, 3, 2, 4, 4, 5, 3, 5, 4, 6, 6, 7, 5, 1, 0, 7, 7, 6, 1, 7, 0, 3, 3, 5, 7, 4, 2, 1, 1, 6, 4])  # A cube
        >>> BSP(igeom=0, ifaces=get_ifaces(0)).invert()
        >>> get_face(0, 0)
        array('i', [2, 1, 0])
        """
        igeom = self.igeom
        # Flip plane
        flip_plane_normal(self.plane)
        # Flip all face normals
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


def clip_ifaces(igeom, ifaces, clipping_bsp, other=None):  # FIXME test
    """
    Recursively remove all ifaces that are inside the clipping_bsp tree.
    """
    master_bsp = geometry[igeom].master_bsp
    plane = clipping_bsp.plane
    front = []
    back = []
    while ifaces:
        iface = ifaces.pop()
        if iface == 49:
            print("here")
        new_coplanar_front, new_coplanar_back, new_front, new_back, spl_ifaces = split_iface(
                igeom=igeom,
                iface=iface,
                plane=plane,  # FIXME move out of for cycle
                )
        front.extend(new_front)
        back.extend(new_back)
        front.extend(new_coplanar_front)
        back.extend(new_coplanar_back)

        # Update lists and bsp tree
        for spl_iface, new_iface in spl_ifaces.items():
            if spl_iface in ifaces:
                ifaces.append(new_iface)
            elif spl_iface in front:
                front.append(new_iface)
            elif spl_iface in back:
                back.append(new_iface)
            elif spl_iface in other:
                print("Found in other!", spl_iface, ">", new_iface)
            else:
                ret = update_iface_in_bsp(master_bsp, spl_iface, new_iface)
                if not ret:
                    print("NOT FOUND iface:", spl_iface, ">", new_iface)
                    if spl_iface == 49:
                        print("here")
#                    raise Exception("Not found in bsp, iface:", igeom, spl_iface)

    if clipping_bsp.front_bsp:
        # recurse on branches, conserve those polygons
        front = clip_ifaces(igeom, front, clipping_bsp.front_bsp, back)

    if clipping_bsp.back_bsp:
        # recurse on branches, conserve those polygons
        back = clip_ifaces(igeom, back, clipping_bsp.back_bsp, front)
    else:  # FIXME this removes clipping capacity
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
            ifaces=bsp.ifaces[:],
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
    return ifaces


def update_iface_in_bsp(bsp, spl_iface, new_iface):  # FIXME test
    """
    Replace iface with new_iface0 and new_iface1 in bsp tree
    """
    if spl_iface in bsp.ifaces:
        bsp.ifaces.append(new_iface)
        return True
    ret_front = False
    ret_back = False
    if bsp.front_bsp:
        ret_front = update_iface_in_bsp(bsp.front_bsp, spl_iface, new_iface)
    if bsp.back_bsp:
        ret_back = update_iface_in_bsp(bsp.back_bsp, spl_iface, new_iface)
    return ret_front or ret_back

# New geom
    

def get_new_geom_from_bsp(bsp):  # FIXME test
#    """
#    Return a new Geom according to bsp contents
#    >>> geometry[0] = Geom([-1.0, -1.0, -1.0, -1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0, 1.0, -1.0, 1.0, 1.0, 1.0, 1.0, 1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0, -1.0], [0, 1, 2, 2, 3, 0, 3, 2, 4, 4, 5, 3, 5, 4, 6, 6, 7, 5, 1, 0, 7, 7, 6, 1, 7, 0, 3, 3, 5, 7, 4, 2, 1, 1, 6, 4])  # A cube
#    >>> bsp = BSP(igeom=0, ifaces=get_ifaces(0))
#    >>> geometry[0] = get_new_geom_from_bsp(bsp)
#    >>> get_nverts(igeom=0), get_nfaces(igeom=0)
#    (8, 12)
#    """
    # Init
    igeom = bsp.igeom
    selected_ifaces = set(get_all_ifaces_from_bsp(bsp))

#    print("len(selected_ifaces) before replacement of parents:",len(selected_ifaces))
#    # Choose the best faces by gruping siblings into parent
#    print("len(selected_ifaces) before:",len(selected_ifaces))
#    done = False
#    while not done:
#        done = True
#        for iface in selected_ifaces:
#            parent = iface_to_parent.get(iface, None)  # could be non existant
#            if parent is None:
#                continue
#            siblings = set(iface_to_children[parent])
#            if siblings and siblings <= selected_ifaces:
#                selected_ifaces -= siblings  # remove siblings
#                selected_ifaces.add(parent)  # and add parent instead
#                done = False  # not finished
#                break  # set is updated, restart of the loop needed
#    print("len(selected_ifaces) after replacement of parents:",len(selected_ifaces))
#
#    print("nverts before:",get_nverts(igeom))
#    merge_duplicated_verts(igeom)
#    print("nverts after:",get_nverts(igeom))
#
#    # Get iverts that deserve a check, those of faces with a parent,
#    # not replaced by the previous step
#    selected_ifaces_w_parent = [iface for iface in selected_ifaces if iface_to_parent.get(iface, None)]
#    check_double_iverts = []
#    for iface in selected_ifaces_w_parent:
#        check_double_iverts.extend(get_face(igeom, iface))
#    # Compare those iverts with the existing ones
#    nverts = get_nverts(igeom)
#    for check_double_ivert in check_double_iverts:
#        # Get the vert coordinates
#        check_double_vert = get_vert(igeom, check_double_ivert)
#        # Compare those with any vert
#        for ivert in range(nverts):
#            if (Vector(get_vert(igeom, ivert))-Vector(check_double_vert)).isZero():
#                # merge, check all faces for that ivert
#                for selected_iface in selected_ifaces:
#                    selected_face = get_face(igeom, selected_iface)
#                    if selected_face[0] == check_double_ivert:
#                        selected_face[0] = ivert
#                        print("Fix ivert:", ivert, check_double_ivert )
#                        break
#                    elif selected_face[1] == check_double_ivert:
#                        selected_face[1] = ivert
#                        print("Fix ivert:", ivert, check_double_ivert )
#                        break
#                    elif selected_face[2] == check_double_ivert:
#                        selected_face[2] = ivert
#                        print("Fix ivert:", ivert, check_double_ivert )
#                        break

# FIXME this break the geometry!!!
#    print("len(selected_ifaces) before merging of triangles:",len(selected_ifaces))
#    # Merge siblings FIXME put in previous routine
#    # The problem are TJunctions.
#    done = False
#    while not done:
#        done = True
#        for iface in selected_ifaces:
#            parent = iface_to_parent.get(iface, None)  # could be non existant
#            if parent is None:  # it is an original face
#                continue
#            descendants = get_iface_descendants(igeom, parent)
#            descendants.remove(iface)
#            for descendant in descendants:
#                if descendant in selected_ifaces:
#                    merged_iface = merge_ifaces(igeom, iface0=iface, iface1=descendant)
##                    print("iface:", iface, "descendant:", descendant, "merged:", merged_iface)
#                    if merged_iface:
#                        # remove merging couple from selected_ifaces
#                        selected_ifaces.remove(iface)
#                        selected_ifaces.remove(descendant)
#                        # insert new iface in selected_ifaces
#                        selected_ifaces.add(merged_iface)
#                        done = False
#                        break
#            if not done:
#                break
#    print("len(selected_ifaces) after merging of triangles:",len(selected_ifaces))

    # Solve local TJunctions FIXME or solve them globally? let's see
    # print("borders:", get_border_halfedges(igeom, selected_ifaces))

#    # Get only selected iverts, avoid duplication
#    selected_iverts = []
#    for iface in selected_ifaces:
#        face = get_face(igeom, iface)
#        selected_iverts.extend((face[0], face[1], face[2]))
#    selected_iverts = set(selected_iverts)
#    # Build the ivert to new_ivert dict, and the corresponding new_verts
#    ivert_to_ivert = {}
#    new_verts = []
#    new_ivert = 0
#    for ivert in selected_iverts:
#        vert = get_vert(igeom, ivert)
#        new_verts.extend(vert)
#        ivert_to_ivert[ivert] = new_ivert
#        new_ivert += 1
#    # Build the new_faces
#    new_faces = []
#    for iface in selected_ifaces:
#        face = get_face(igeom, iface)
#        new_faces.extend((
#                ivert_to_ivert[face[0]],
#                ivert_to_ivert[face[1]],
#                ivert_to_ivert[face[2]],
#                ))
#    print("nverts after removal of duplicates:", len(new_verts)/3)
    
    # The simplest builder
    new_verts = geometry[igeom].verts
    new_faces = []
    for selected_iface in selected_ifaces:
        face = get_face(igeom, selected_iface)
        new_faces.extend(face)
    return Geom(new_verts, new_faces)


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
#    clip_to(b, a)  # remove everything in b inside a
#    b.invert()
#    clip_to(b, a)  # remove everything in -b inside a
#    b.invert()

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


if __name__ == "__main__":
    import doctest
    doctest.testmod()