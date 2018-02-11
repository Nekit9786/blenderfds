#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 10 17:12:27 2018

@author: egissi
"""
import array
from vector import Vector

EPSILON = 1e-07
EPSILON2 = 1e-05

geometry = [None, None, None, None]


# Notation

# geometry: list of Geom()
# igeom: 2, index of a Geom() in geometry

# faces: [0,1,2, 0,2,3], list of connectivity of all faces
# nfaces: range(int(len(faces)/3))
# face: [0,2,3], connectivity of one face

# ifaces: [0,1,2,7,11,], list of indexes of selected faces
# iface: 7, index of a face

# verts: [-1,-1,1, 1,-1,1, 1,1,1, -1,1,1], list of coo of all vertices
# nverts: range(int(len(verts)/3))
# vert: [1,-1,1], Vector of one vertex

# iverts: [3,4,5,9,], list of indexes of selected verts
# ivert: 7, index of a vert

# isurfids: [0,1,1,1,4], list of isurfid for each face
# isurfid: 4, index of SURF_IDV


class Geom():
    def __init__(self, verts, faces, isurfids=None):
        # Check lenght
        if (len(verts) % 3) != 0:
            raise Exception('Invalid geom, verts lenght not 3n.')
        if (len(faces) % 3) != 0:
            raise Exception('Invalid geom, faces lenght not 3n.')
        # Init tmp geometry
        self.verts = array.array('f', verts)
        self.faces = array.array('i', faces)
        # Keep a link between ifaces and their parent
        # Example: {3: 1, 4: 1} ifaces 3 and 4 are children of 1
        self.iface_to_parent = {}
        # Keep a reference of the original face surf_idi, surf_id index
        # Not updated, only original faces have it
        if isurfids:
            if len(isurfids) != (len(faces) / 3):
                raise Exception('Invalid geom, isurfids lenght not corresponding to nfaces.')
            else:
                self.isurfids = array.array('i', isurfids)
        else:  # Same as iface
            self.isurfids = array.array('i', range((len(faces) // 3)))

#    def __repr__(self):
#        return 'Geom(\n     {},\n     {},\n)'.format(self.verts, self.faces,)

    def __repr__(self):
        text = 'Geom:'
        for iface in range(int(len(self.faces)/3)):
            face = self.faces[3*iface:3*iface+3]
            text += '\n{0}-{4}: {1[0]:.1f},{1[1]:.1f},{1[2]:.1f}, {2[0]:.1f},{2[1]:.1f},{2[2]:.1f}, {3[0]:.1f},{3[1]:.1f},{3[2]:.1f}'.format(
                    iface,
                    self.get_vert(face[0]),
                    self.get_vert(face[1]),
                    self.get_vert(face[2]),
                    self.get_isurfid(iface),
                    )
        return text

    def clone(self):
        return Geom(verts=self.verts[:], faces=self.faces[:])

    def get_vert(self, ivert):  # FIXME duplicated in functional!
        return self.verts[3*ivert:3*ivert+3]

    def get_isurfid(self, iface):  # FIXME duplicated in functional!
        try:
            return self.isurfids[iface]
        except IndexError:
            return None


# iface parent, children


def get_iface_to_parent(igeom):
    """
    Get iface to parent dict
    >>> geometry[0] = Geom([], [0,1,3, 0,1,4, 1,3,4, 1,3,4, 1,3,4, 1,3,4, 1,3,4, 1,3,4, 1,3,4,])
    >>> geometry[0].iface_to_parent = {2:0, 3:2, 4:2, 5:1, 6:3, 7:6, 8:3}
    >>> get_iface_to_parent(igeom=0)
    {2: 0, 3: 2, 4: 2, 5: 1, 6: 3, 7: 6, 8: 3}
    """
    return geometry[igeom].iface_to_parent


def get_iface_to_children(igeom):
    """
    Build iface_to_children dict
    >>> geometry[0] = Geom([], [0,1,3, 0,1,4, 1,3,4, 1,3,4, 1,3,4, 1,3,4, 1,3,4, 1,3,4, 1,3,4,])
    >>> geometry[0].iface_to_parent = {2:0, 3:2, 4:2, 5:1, 6:3, 7:6, 8:3}
    >>> get_iface_to_children(igeom=0)
    {0: [2], 2: [3, 4], 1: [5], 3: [6, 8], 6: [7]}
    """
    iface_to_children = {}
    iface_to_parent = get_iface_to_parent(igeom)
    for child, parent in iface_to_parent.items():
        try:
            iface_to_children[parent].append(child)
        except KeyError:
            iface_to_children[parent] = [child, ]
    return iface_to_children


def get_iface_descendants(igeom, iface, iface_to_children=None):
    """
    Return all iface descendants
    >>> geometry[0] = Geom([], [0,1,3, 0,1,4, 1,3,4, 1,3,4, 1,3,4, 1,3,4, 1,3,4, 1,3,4, 1,3,4,])
    >>> geometry[0].iface_to_parent = {2:0, 3:2, 4:2, 5:1, 6:3, 7:6, 8:3}
    >>> get_iface_descendants(igeom=0, iface=0)
    [2, 3, 4, 6, 8, 7]
    """
    # Init
    if iface_to_children is None:
        iface_to_children = get_iface_to_children(igeom)
    # Children
    children = iface_to_children.get(iface, [])
    descendants = children[:]
    for child in children:
        descendants.extend(get_iface_descendants(
                igeom, child, iface_to_children
                ))
    return descendants


# faces


def get_face(igeom, iface):
    """
    Get iface face connectivity
    >>> geometry[0] = Geom([-1,-1,1, 1,-1,1, 1,1,1, -1,1,1], [0,1,2, 0,2,3])
    >>> get_face(0, 1)
    array('i', [0, 2, 3])
    """
    return geometry[igeom].faces[3*iface:3*iface+3]


def set_face(igeom, iface, face):
    """
    Set iface face connectivity
    >>> geometry[0] = Geom([-1,-1,1, 1,-1,1, 1,1,1, -1,1,1], [0,1,2, 0,2,3])
    >>> set_face(0, 1, (3,2,0))
    >>> get_face(0, 1)
    array('i', [3, 2, 0])
    """
    geometry[igeom].faces[3*iface:3*iface+3] = array.array('i', face)


def append_face(igeom, face, iface_parent=None):
    """
    Append a face to the Geom, register its parent, return its iface.
    >>> geometry[0] = Geom([-1,-1,1, 1,-1,1, 1,1,1, -1,1,1], [0,1,2, ])
    >>> iface = append_face(0, [0,2,3], 0)
    >>> get_face(0, iface)
    array('i', [0, 2, 3])
    """
    # Append face
    geometry[igeom].faces.extend(face)
    iface = get_nfaces(igeom) - 1
    # Register parent
    if iface_parent is not None:
        geometry[igeom].iface_to_parent[iface] = iface_parent
    # Return
    return iface


def get_nfaces(igeom):
    """
    Get the len of faces
    >>> geometry[0] = Geom([-1,-1,1, 1,-1,1, 1,1,1, -1,1,1], [0,1,2, 0,2,3])
    >>> get_nfaces(0)
    2
    """
    return int(len(geometry[igeom].faces)/3)


def get_ifaces(igeom):
    """
    Get the range for faces
    >>> geometry[0] = Geom([-1,-1,1, 1,-1,1, 1,1,1, -1,1,1], [0,1,2, 0,2,3])
    >>> get_ifaces(0)
    [0, 1]
    """
    return [i for i in range(int(len(geometry[igeom].faces)/3))]


def get_iface_plane(igeom, iface):
    """
    Get plane from iface normal and distance.
    >>> geometry[0] = Geom([1,0,0, 0,1,0, 0,0,1,], [0,1,2, ])
    >>> n, w = get_iface_plane(0, 0)
    >>> print('n:', n, 'w:', w)
    n: Vector(0.577, 0.577, 0.577) w: 0.5773502691896258
    """
    face = get_face(igeom, iface)
    a, b, c = (
        Vector(get_vert(igeom, face[0])),
        Vector(get_vert(igeom, face[1])),
        Vector(get_vert(igeom, face[2])),
    )
    normal = b.minus(a).cross(c.minus(a)).unit()  # plane normal vector
    distance = normal.dot(a)  # plane distance from origin
    return normal, distance


def flip_iface_normal(igeom, iface):
    """
    Flip face normal.
    >>> geometry[0] = Geom([1,0,0, 0,1,0, 0,0,1,], [0,1,2, ])
    >>> flip_iface_normal(igeom=0, iface=0)
    >>> n, w = get_iface_plane(igeom=0, iface=0)
    >>> print('n:', n, 'w:', w)
    n: Vector(-0.577, -0.577, -0.577) w: -0.5773502691896258
    """
    face = get_face(igeom, iface)
    set_face(igeom, iface, face=(face[2], face[1], face[0]))


def split_iface(igeom, iface, spl_igeom, spl_iface):
    """
    Split iface from igeom by spl_iface of spl_igeom if needed.
    Append new verts and new faces to geometry igeom.
    Return ifaces in the appropriate lists.
    >>> geometry[0] = Geom([0.0,0.0,1.0, 0.0,0.0,-1.0, 0.0,1.0,1.0],[0,1,2]) # axis +x
    >>> geometry[1] = Geom([-1.0,0.0,0.0, 1.0,0.0,0.0, 0.0,1.0,0.0],[0,1,2]) # axis +z
    >>> geometry[2] = Geom([-1.0,0.0,0.0, 1.0,0.0,0.0, 0.0,1.0,0.0],[2,1,0]) # axis -z
    >>> split_iface(igeom=0, iface=0, spl_igeom=1, spl_iface=0)
    ([], [], [1, 2], [3], [3, 4])
    >>> geometry[0]
    Geom:
    0-0: 0.0,0.0,1.0, 0.0,0.0,-1.0, 0.0,1.0,1.0
    1-None: 0.0,0.0,1.0, 0.0,0.0,0.0, 0.0,0.5,0.0
    2-None: 0.0,0.0,1.0, 0.0,0.5,0.0, 0.0,1.0,1.0
    3-None: 0.0,0.0,0.0, 0.0,0.0,-1.0, 0.0,0.5,0.0
    >>> split_iface(igeom=1, iface=0, spl_igeom=0, spl_iface=2)
    ([], [], [1], [2], [3])
    >>> geometry[1]
    Geom:
    0-0: -1.0,0.0,0.0, 1.0,0.0,0.0, 0.0,1.0,0.0
    1-None: 0.0,0.0,0.0, 1.0,0.0,0.0, 0.0,1.0,0.0
    2-None: -1.0,0.0,0.0, 0.0,0.0,0.0, 0.0,1.0,0.0
    >>> split_iface(igeom=1, iface=1, spl_igeom=1, spl_iface=1)  # coplanar front
    ([1], [], [], [], [])
    >>> split_iface(igeom=1, iface=1, spl_igeom=2, spl_iface=0)  # coplanar back
    ([], [1], [], [], [])
    """
    # Vertices and faces types, collections of ifaces
    COPLANAR = 0  # vertex or face is within EPSILON2 distance from plane
    FRONT = 1     # vertex or face is in front of the plane
    BACK = 2      # vertex or face is at the back of the plane
    SPANNING = 3  # edge is intersected
    coplanar_front, coplanar_back, front, back = [], [], [], []
    cut_iverts = []  # collection of new iverts due to cuts

    # Classify each point as well as the entire polygon
    # into one of the above classes.
    faceType = 0
    vertexTypes = []
    spl_normal, spl_distance = get_iface_plane(spl_igeom, spl_iface)
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
        iface_normal, iface_w = get_iface_plane(igeom, iface)
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
                # Record cut_ivert for later fixing of TJunctions FIXME
                cut_iverts.append(cut_ivert)

        # Add front cut faces
        if len(front_iverts) > 2:
            new_iface = append_face(igeom, front_iverts[0:3], iface)
            front.append(new_iface)
            if len(front_iverts) == 4:
                new_iface = append_face(
                        igeom,
                        (front_iverts[0], front_iverts[2], front_iverts[3]),
                        iface,
                    )
                front.append(new_iface)
        else:
            raise Exception('Problem with front spanning:', front_iverts)

        # Add back cut faces
        if len(back_iverts) > 2:
            new_iface = append_face(igeom, back_iverts[0:3], iface)
            back.append(new_iface)
            if len(back_iverts) == 4:
                new_iface = append_face(
                        igeom,
                        (back_iverts[0], back_iverts[2], back_iverts[3]),
                        iface,
                    )
                back.append(new_iface)
        else:
            raise Exception('Problem with back spanning:', back_iverts)

    # Return
    return coplanar_front, coplanar_back, front, back, cut_iverts


def merge_ifaces(igeom, iface0, iface1):  # FIXME broken normals?
    """
    Merge two ifaces from igeom
    Return new iface
    >>> geometry[0] = Geom([-1,0,0, 0,0,0, 1,0,0, 0,1,0], [0,1,3, 1,2,3])  # 012 aligned, +z
    >>> iface = merge_ifaces(igeom=0, iface0=0, iface1=1)
    >>> get_face(igeom=0, iface=iface)
    array('i', [0, 2, 3])
    >>> geometry[0] = Geom([-1,0,0, 0,0,0, 1,0,0, 0,1,0], [0,3,1, 1,3,2])  # 012 aligned, -z
    >>> iface = merge_ifaces(igeom=0, iface0=0, iface1=1)
    >>> get_face(igeom=0, iface=iface)
    array('i', [0, 3, 2])
    >>> geometry[0] = Geom([-1,0,0, 0,-1,0, 1,0,0, 0,0,0], [0,1,3, 1,2,3])  # 032 aligned, +z
    >>> iface = merge_ifaces(igeom=0, iface0=0, iface1=1)
    >>> get_face(igeom=0, iface=iface)
    array('i', [0, 1, 2])
    >>> geometry[0] = Geom([-1,0,0, 0,-1,0, 1,0,0, 0,0,0], [0,3,1, 1,3,2])  # 032 aligned, -z
    >>> iface = merge_ifaces(igeom=0, iface0=0, iface1=1)
    >>> get_face(igeom=0, iface=iface)
    array('i', [0, 2, 1])
    >>> geometry[0] = Geom([-1,0,0, 0,-1,0, 1,0,0, 0,0,0], [3,1,0, 1,3,2])  # 032 aligned, -z, reorder
    >>> iface = merge_ifaces(igeom=0, iface0=0, iface1=1)
    >>> get_face(igeom=0, iface=iface)
    array('i', [0, 2, 1])
    >>> geometry[0] = Geom([-1,0,0, 0,-1,0, 1,0,0, 0,0,0], [3,1,0, 2,1,3])  # 032 aligned, -z, reorder
    >>> iface = merge_ifaces(igeom=0, iface0=0, iface1=1)
    >>> get_face(igeom=0, iface=iface)
    array('i', [0, 2, 1])
    """
    # Get parent
    iface_to_parent = get_iface_to_parent(igeom)
    iface_parent = iface_to_parent.get(iface0, None)
    # calc new face and append it to geometry[igeom]
    # update geometry[igeom].iface_to_parent
    # Get shared and unshared vertices
    face0 = get_face(igeom, iface0)
    face1 = get_face(igeom, iface1)
    shared_iverts = []
    unshared_iverts = []
    for ivert in face0:
        if ivert in face1:
            shared_iverts.append(ivert)
        else:
            unshared_iverts.append(ivert)
    for ivert in face1:
        if ivert not in face0:
            unshared_iverts.append(ivert)
    # Merge if possible
    if len(shared_iverts) == 2:
        v0 = get_vert(igeom, unshared_iverts[0])
        v1 = get_vert(igeom, unshared_iverts[1])
        for i, shared_ivert in enumerate(shared_iverts):
            vu = get_vert(igeom, shared_ivert)
            if vu.is_within(v0, v1) and vu.is_collinear(v0, v1):
                if i == 0:
                    face = (
                            unshared_iverts[0],
                            unshared_iverts[1],
                            shared_iverts[1],
                            )
                else:
                    face = (
                            unshared_iverts[0],
                            shared_iverts[0],
                            unshared_iverts[1],
                            )
                iface = append_face(igeom, face, iface_parent)
                iface_to_parent[iface0] = iface  # FIXME
                iface_to_parent[iface1] = iface  # FIXME
                return iface


# verts


def get_vert(igeom, ivert):
    """
    Get ivert vertex coordinates
    >>> geometry[0] = Geom([-1,-1,1, 1,-1,1, 1,1,1, -1,1,1], [0,1,2, 0,2,3])
    >>> get_vert(0, 2)
    Vector(1.000, 1.000, 1.000)
    """
    return Vector(geometry[igeom].verts[3*ivert:3*ivert+3])


def append_vert(igeom, vert):
    """
    Append a vert to the Geom, return its index.
    >>> geometry[0] = Geom([-1,-1,1, 1,-1,1, 1,1,1, -1,1,1], [0,1,2, 0,2,3])
    >>> ivert = append_vert(0, Vector(1., 1., 2.))
    >>> get_vert(0, ivert)
    Vector(1.000, 1.000, 2.000)
    """
    geometry[igeom].verts.extend(list(vert))
    return get_nverts(igeom)-1


def get_nverts(igeom):
    """
    Get the len of vertices
    >>> geometry[0] = Geom([-1,-1,1, 1,-1,1, 1,1,1, -1,1,1], [0,1,2, 0,2,3])
    >>> get_nverts(0)
    4
    """
    return int(len(geometry[igeom].verts)/3)


def get_iverts(igeom):
    """
    Get the range for verts
    >>> geometry[0] = Geom([-1,-1,1, 1,-1,1, 1,1,1, -1,1,1], [0,1,2, 0,2,3])
    >>> get_iverts(0)
    [0, 1, 2, 3]
    """
    return [i for i in range(int(len(geometry[igeom].verts)/3))]


def merge_duplicated_verts(igeom):
    """
    Remove dup verts, and relink all faces. No mod of faces.
    >>> geometry[0] = Geom([-1,0,0, 0,0,0, 1,0,0, 0,1,0, 0,1,0], [0,1,3, 0,1,4, 1,3,4,])
    >>> merge_duplicated_verts(igeom=0)
    >>> geometry[0].faces
    array('i', [0, 1, 3, 0, 1, 3, 1, 3, 3])
    >>> geometry[0].verts
    array('f', [-1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0])
    """
    # Find unique verts and build ivert_to_ivert dict
    selected_pyverts = []
    ivert_to_ivert = {}
    nverts = get_nverts(igeom)
    for ivert in range(nverts):
        vert = Vector(get_vert(igeom, ivert))
        seen = False
        for i, selected_pyvert in enumerate(selected_pyverts):
            if (selected_pyvert - vert).isZero():
                seen = True
                ivert_to_ivert[ivert] = i
                break
        if not seen:
            selected_pyverts.append(vert)
            ivert_to_ivert[ivert] = len(selected_pyverts) - 1
    # Update face ivert links
    for i, ivert in enumerate(geometry[igeom].faces):
        new_ivert = ivert_to_ivert[ivert]
        geometry[igeom].faces[i] = new_ivert
    # Flat pyverts to build new verts
    selected_verts = array.array('f')
    for pyvert in selected_pyverts:
        selected_verts.extend(pyvert)
    geometry[igeom].verts = selected_verts


# edges


def get_halfedges(igeom, ifaces):
    """
    Get halfedges dict
    Eg: {(1,2):7]} with {(ivert0, ivert1) : iface on the left}
    according to iface0 normal up
    >>> geometry[0] = Geom([-1,-1,1, 1,-1,1, 1,1,1, -1,1,1], [0,1,2, 0,2,3])  # open
    >>> print(get_halfedges(igeom=0, ifaces=get_ifaces(igeom=0)))
    {(0, 1): 0, (1, 2): 0, (2, 0): 0, (0, 2): 1, (2, 3): 1, (3, 0): 1}
    >>> geometry[0] = Geom([-1,-1,1, 1,-1,1, 1,1,1, -1,1,1], [0,1,2, 0,3,2])  # open, unorientable
    >>> print(get_halfedges(igeom=0, ifaces=get_ifaces(igeom=0)))
    Traceback (most recent call last):
    ...
    Exception: ('Invalid GEOM, non-manifold or unorientable at iface:', 1)
    """
    halfedges = dict()
    for iface in ifaces:
        face = get_face(igeom, iface)
        for i in range(3):
            halfedge = (face[i], face[(i+1) % 3])
            if halfedge in halfedges:
                raise Exception('Invalid GEOM, non-manifold or unorientable at iface:', iface)
            halfedges[halfedge] = iface
    return halfedges


def get_border_halfedges(igeom, ifaces):
    """
    Get border halfedges dict
    Eg: {(1,2):7]} with {(ivert0, ivert1): iface on the left}
    according to iface0 normal up
    >>> geometry[0] = Geom([-1,-1,1, 1,-1,1, 1,1,1, -1,1,1], [0,1,2, 0,2,3])  # open
    >>> print(get_border_halfedges(igeom=0, ifaces=get_ifaces(igeom=0)))
    {(3, 0): 1, (2, 3): 1, (1, 2): 0, (0, 1): 0}
    """
    halfedges = get_halfedges(igeom, ifaces)
    border_halfedges = {}
    while halfedges:
        halfedge, iface = halfedges.popitem()
        opposite = halfedge[1], halfedge[0]
        if opposite in halfedges:
            del halfedges[opposite]
        else:
            border_halfedges[halfedge] = iface
    return border_halfedges


def get_edge_loops(igeom, border_halfedges):  # FIXME working here
    """
    Get oriented edge loops,
    Eg: [3,0,1,2,] with ivert0, ivert1, ...
    according to iface0 normal up
    >>> geometry[0] = Geom([-1,-1,1, 1,-1,1, 1,1,1, -1,1,1], [0,1,2, 0,2,3])  # open
    >>> border_halfedges = get_border_halfedges(igeom=0, ifaces=get_ifaces(igeom=0))
    >>> get_edge_loops(0, border_halfedges)
    [[0, 1, 2, 3]]
    >>> geometry[0] = Geom([-1,0,0, 0,0,0, 1,0,0, 1,1,0, -1,1,0, -1,-1,0, 0,-1,0, 1,-1,0],\
                [1,3,4, 1,2,3, 0,5,6, 6,2,0, 6,7,2, 0,1,4,])  # open
    >>> border_halfedges = get_border_halfedges(igeom=0, ifaces=get_ifaces(igeom=0))
    >>> get_edge_loops(0, border_halfedges)
    [[3, 4, 0, 5, 6, 7, 2], [1, 2, 0]]
    """
    border_halfedges = list(border_halfedges)  # FIXME
    edge_loops = []
    while border_halfedges:
        bh = border_halfedges.pop()
        edge_loop = [bh[0], bh[1]]
        while edge_loop[-1] != edge_loop[0]:
            # Get candidates
            candidates = []
            for bh in border_halfedges:
                if edge_loop[-1] == bh[0]:
                    candidates.append(bh)
            # Good or singularity? Choose the candidate that goes back home
            if len(candidates) == 1:
                choosen = 0
            elif len(candidates) == 2:
                ve = Vector(get_vert(igeom, edge_loop[-2])) \
                    - Vector(get_vert(igeom, edge_loop[-1]))
                dot0 = Vector(get_vert(igeom, candidates[0][0])).dot(ve)
                dot1 = Vector(get_vert(igeom, candidates[1][0])).dot(ve)
                if dot0 < dot1:
                    choosen = 0
                else:
                    choosen = 1
            else:  # In case of 0 or more than 2 candidates
                raise Exception('Invalid GEOM, non-manifold in edge_loop:', edge_loop)
            # Append candidate
            border_halfedges.remove(candidates[choosen])
            edge_loop.append(candidates[choosen][1])
        # Closed
        edge_loop.pop()  # Remove last duplicated
        edge_loops.append(edge_loop[:])
    return edge_loops


# STL


def to_STL(igeom, filename):
    """
    Write self to STL file
    """
    with open(filename, 'w') as f:
        f.write('solid name\n')
        for iface in get_ifaces(igeom):
            f.write('facet normal 0 0 0\n')
            f.write(' outer loop\n')
            for ivert in get_face(igeom, iface):
                f.write('  vertex {v[0]:.6f} {v[1]:.6f} {v[2]:.6f}\n'.format(
                        v=get_vert(igeom, ivert)))
            f.write(' endloop\n')
            f.write('endfacet\n')
        f.write('endsolid name\n')
    print('to_STL:', filename)


def from_STL(filename):
    """
    Import verts and faces from STL file into self
    >>> geometry[0] = Geom([-1,-1,1, 1,-1,1, 1,1,1, -1,1,1], [0,1,2, 0,2,3])
    >>> to_STL(0,'./test/doctest.stl')
    to_STL: ./test/doctest.stl
    >>> from_STL('./test/doctest.stl')
    Geom:
    0-0: -1.0,-1.0,1.0, 1.0,-1.0,1.0, 1.0,1.0,1.0
    1-1: -1.0,-1.0,1.0, 1.0,1.0,1.0, -1.0,1.0,1.0
    """
    # Get STL mesh
    from stl import mesh
    mesh = mesh.Mesh.from_file(filename)
    verts, faces, py_verts = [], [], []
    for iface, p in enumerate(mesh.points):
        # p is [-1.0, -1.0, -1.0, -1.0, -1.0, 1.0, -1.0, 1.0, 1.0]
        verts.extend(p)
        faces.extend((3*iface, 3*iface+1, 3*iface+2))
        py_verts.append((p[0], p[1], p[2]))
        py_verts.append((p[3], p[4], p[5]))
        py_verts.append((p[6], p[7], p[8]))
    # Remove duplicated verts
    unique_py_verts = [x for n, x in enumerate(py_verts)
                       if x not in py_verts[:n]]
    new_faces = []
    for ivert in faces:
        new_ivert = unique_py_verts.index(py_verts[ivert])  # Get new index
        new_faces.append(new_ivert)  # Append it to the faces
    unique_verts = [v for vs in unique_py_verts for v in vs]  # Flatten
    return Geom(unique_verts, new_faces)


# Geometry sanity

def check_loose_verts(igeom):
    """
    Check loose vertices: vertices that have no connectivity
    >>> geometry[0] = Geom([-1,-1,0, 1,-1,0, 0,1,0, 0,0,1,  0,0,2], [2,1,0, 0,1,3, 1,2,3, 2,0,3])  # Tet, loose vert
    >>> check_loose_verts(igeom=0)
    Traceback (most recent call last):
    ...
    Exception: Invalid GEOM, loose verts.
    """
    nverts = get_nverts(igeom)
    nfaces = get_nfaces(igeom)
    used_iverts = []
    for iface in range(nfaces):
        used_iverts.extend(get_face(igeom, iface))
    used_iverts = set(used_iverts)
    if nverts != len(used_iverts) or nverts != max(used_iverts) + 1:
        raise Exception("Invalid GEOM, loose verts.")


def check_degenerate_geometry(igeom):
    """
    Check degenerate geometry: zero lenght edge, zero area faces
    >>> geometry[0] = Geom([-1,-1,0, 1,-1,0, 0,0,1, 0,0,1], [2,1,0, 0,1,3, 1,2,3, 2,0,3])  # Tet, zero edge
    >>> check_degenerate_geometry(igeom=0)
    Traceback (most recent call last):
    ...
    Exception: ('Invalid GEOM, zero lenght edge in face:', 2)
    >>> geometry[0] = Geom([-1,-1,0, 1,-1,0, 0,-1,0, 0,0,1], [2,1,0, 0,1,3, 1,2,3, 2,0,3])  # Tet, zero face
    >>> check_geom_sanity(0)
    Traceback (most recent call last):
    ...
    Exception: ('Invalid GEOM, zero area iface:', 0)
    """
    nfaces = get_nfaces(igeom)
    for iface in range(nfaces):
        face = get_face(igeom, iface)
        a, b, c = \
            get_vert(igeom, face[0]), \
            get_vert(igeom, face[1]), \
            get_vert(igeom, face[2])
        if (a - b).isZero() or (b - c).isZero() or (c - a).isZero():
            raise Exception("Invalid GEOM, zero lenght edge in face:", iface)
        if (a - b).cross(a - c).isZero():
            raise Exception('Invalid GEOM, zero area iface:', iface)


def check_is_solid(igeom):
    """
    Check surface:
    - 2-manifold and closed, each edge should join two faces, no more no less
    - orientable, adjoining faces should have normals in the same directions
    >>> geometry[0] = Geom([-1,-1,0, 1,-1,0, 0,1,0, 0,0,1], [       0,1,3, 1,2,3, 2,0,3])  # Tet, no base
    >>> check_is_solid(0)
    Traceback (most recent call last):
    ...
    Exception: ('Invalid GEOM, non closed at ifaces:', [2, 1, 0])
    """
    border_halfedges = get_border_halfedges(igeom, get_ifaces(igeom))
    if border_halfedges:
        raise Exception("Invalid GEOM, non closed at ifaces:",
                        [b for b in border_halfedges.values()])

def check_euler(igeom):  # FIXME test
    """
    Euler formula: nverts - nedges + nfaces = chi
    Euler characteristic chi of the connected sum of g tori is:
    chi = 2 − 2g, with g genus
    g = 0, 1, 2, 3, ... => chi = 2, 0, -2, -4, ...
    """
    nverts = get_nverts(igeom)
    nfaces = get_nfaces(igeom)
    # Count edges
    nedges = 0
    seen_halfedges = []
    halfedges = get_halfedges(igeom, get_ifaces(igeom))
    for halfedge in halfedges:
        opposite = halfedge[1], halfedge[0]
        if halfedge in seen_halfedges or opposite in seen_halfedges:
            continue
        else:
            seen_halfedges.append(halfedge)
            nedges += 1
    # Check Euler formula
    chi = nverts - nedges + nfaces
    if chi not in range(2, 100, 2):
        raise Exception('Invalid GEOM, chi in Euler formula is:', chi)


def check_geom_sanity(igeom):
    """
    Check geometry sanity

    If the mesh is correct and encloses a volume, this can be checked with
    prior tests: checking orientability, non-borders, non-self-intersecting.
    After that we can calculate its topological features and check if
    Euler's formula  C+V=A+2(S-H) is satisfied.
    If the mesh is not correct, many geometric algorithms will fail.
    The only solution in this case is the user repairing the mesh.
    >>> geometry[0] = Geom([-1,-1,0, 1,-1,0, 0,1,0, 0,0,1], [2,1,0, 0,1,3, 1,2,3, 2,0,3])  # Good tet
    >>> check_geom_sanity(0)
    True

#    >>> geometry[0] = from_STL('self_intersecting.stl')
#    >>> check_geom_sanity(0)
#    Traceback (most recent call last):
#    ...
#    Exception: ('Invalid GEOM, self-intersecting:', 79, 145)
    """
    check_loose_verts(igeom)
    check_degenerate_geometry(igeom)
    check_is_solid(igeom)
    check_euler(igeom)

    # Check correct normals for a solid in fluid FIXME working here

    # Check self intersection  # FIXME not working
#    ifaces = get_ifaces(igeom)
#    bsp = build_bsp(igeom, ifaces)
#    inverted_bsp = get_inverted_bsp(bsp)
#    clip_to(bsp, inverted_bsp)
#    remaining_ifaces = get_all_ifaces_from_bsp(bsp)
#    if remaining_ifaces:
#        raise Exception('Invalid GEOM, self-intersecting at faces:', remaining_ifaces)

    return True


if __name__ == "__main__":
    import doctest
    doctest.testmod()
