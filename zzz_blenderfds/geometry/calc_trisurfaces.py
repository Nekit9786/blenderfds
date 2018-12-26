"""BlenderFDS, algorithms for triangulated surfaces."""

import bpy, bmesh, mathutils
from time import time
from math import floor, ceil

from ..exceptions import BFException
from . import utils

DEBUG = False


# Get triangulated surface

def get_trisurface(context, ob, check=True) -> "mas, verts, faces":
    """Get triangulated surface from object ready for FDS GEOM format."""
    # Check and init
    DEBUG and print("BFDS: get_triangles")
    assert(ob.type == 'MESH')
    if not ob.data.vertices:
        raise BFException(ob, "Empty object!")
    # Check original mesh quality
    if check:
        check_mesh_quality(context, ob)
    # Create new object global copy
    ob_tmp = utils.object_get_global_copy(context, ob, suffix='_tri_tmp')
    # Create triangulate modifier
    mo = ob_tmp.modifiers.new('triangulate_tmp','TRIANGULATE')
    mo.quad_method, mo.ngon_method = 'BEAUTY', 'BEAUTY'
    # Apply triangulate modifier and remove it
    ob_tmp.data = ob_tmp.to_mesh(
        scene=context.scene,
        apply_modifiers=True,
        settings="RENDER",
        calc_tessface=True,
        calc_undeformed=False,
    )
    ob_tmp.modifiers.remove(mo)
    # Get ob materials from slots
    mas = list()
    material_slots = ob.material_slots
    if len(material_slots) == 0:
        bpy.data.objects.remove(ob_tmp, True)
        raise BFException(ob,
            "No referenced SURF, add at least one Material")
    for material_slot in material_slots:
        ma = material_slot.material
        if not ma:
            bpy.data.objects.remove(ob_tmp, True)
            raise BFException(ob,
                "No referenced SURF, fill empty slot with Material")
        if not ma.bf_export:
            bpy.data.objects.remove(ob_tmp, True)
            raise BFException(ob,
                "Referenced SURF ID='{}' is not exported".format(ma.name))
        mas.append(ma.name)
    # Get ob verts and faces
    bm = bmesh.new()
    bm.from_mesh(ob_tmp.data)
    verts = [(v.co.x, v.co.y, v.co.z) for v in bm.verts]
    faces = [(
            f.verts[0].index+1, # FDS index start from 1, not 0
            f.verts[1].index+1,
            f.verts[2].index+1,
            f.material_index+1,
    ) for f in bm.faces]
    # Clean up
    bm.free()
    bpy.data.objects.remove(ob_tmp, True)
    return mas, verts, faces


# Check mesh quality

def check_mesh_quality(context, ob):
    """Check that Object is a closed orientable manifold,
    with no degenerate geometry."""
    # Init
    DEBUG and print("BFDS: check_mesh_quality")
    bpy.ops.object.mode_set(mode='OBJECT')
    bm = bmesh.new()
    bm.from_mesh(ob.data)
    bm.faces.ensure_lookup_table()  # update bmesh index
    bm.verts.ensure_lookup_table()  # update bmesh index
    bm.edges.ensure_lookup_table()  # update bmesh index
    epsilon_len = context.scene.bf_config_min_edge_length
    epsilon_area = context.scene.bf_config_min_face_area
    bad_verts, bad_edges, bad_faces= list(), list(), list()
    # Check manifold edges, each edge should join two faces, no more no less
    for edge in bm.edges:
        if not edge.is_manifold:
            bad_edges.append(edge)
    if bad_edges:
        msg = "Non manifold or open geometry detected, bad edges selected."
        _raise_bad_geometry(context, ob, bm, msg, bad_edges=bad_edges)
    # Check manifold vertices
    for vert in bm.verts:
        if not vert.is_manifold:
            bad_verts.append(vert)
    if bad_verts:
        msg = "Non manifold vertices detected, bad vertices selected."
        _raise_bad_geometry(context, ob, bm, msg, bad_verts=bad_verts)
    # Check contiguous normals, adjoining faces should have normals
    # in the same directions
    for edge in bm.edges:
        if not edge.is_contiguous:
            bad_edges.append(edge)
    if bad_edges:
        msg = "Inconsistent face normals detected, bad edges selected."
        _raise_bad_geometry(context, ob, bm, msg, bad_edges=bad_edges)
    # Check no degenerate edges, zero lenght edges
    for edge in bm.edges:
        if edge.calc_length() <= epsilon_len:
            bad_edges.append(edge)
    if bad_edges:
        msg = "Too short edges detected, bad edges selected."
        _raise_bad_geometry(context, ob, bm, msg, bad_edges=bad_edges)
    # Check degenerate faces, zero area faces
    for face in bm.faces:
        if face.calc_area() <= epsilon_area:
            bad_faces.append(face)
    if bad_faces:
        msg = "Too small area faces detected, bad faces selected."
        _raise_bad_geometry(context, ob, bm, msg, bad_faces=bad_faces)
    # Check loose vertices, vertices that have no connectivity
    for vert in bm.verts:
        if not bool(vert.link_edges):
            bad_verts.append(vert)
    if bad_verts:
        msg = "Loose vertices detected, bad vertices selected."
        _raise_bad_geometry(context, ob, bm, msg, bad_verts=bad_verts)
    # Check duplicate vertices
    _check_duplicate_vertices(context, ob, bm, epsilon_len)
    # Set overall normals
    bmesh.ops.recalc_face_normals(bm, faces=bm.faces)
    bm.to_mesh(ob.data)
    # Close
    bm.free()

def _check_duplicate_vertices(context, ob, bm, epsilon_len):
    """Check duplicate vertices."""
    bad_verts = list()
    size = len(bm.verts)
    kd = mathutils.kdtree.KDTree(size)  # create a kd-tree from a mesh
    for i, vert in enumerate(bm.verts):
        kd.insert(vert.co, i)
    kd.balance()
    for vert in bm.verts:
        vert_group = list()
        for (co, i, dist) in kd.find_range(vert.co, epsilon_len):
            vert_group.append(i)
        if len(vert_group) > 1:
            for i in vert_group:
                bad_verts.append(bm.verts[i])
    if bad_verts:
        msg = "Duplicate vertices detected, bad vertices selected."
        _raise_bad_geometry(context, ob, bm, msg, bad_verts=bad_verts)


# Check intersections

def _get_bm_and_tree(context, ob, epsilon_len=0., matrix=None):
    """Get BMesh and BVHTree from Object."""
    bm = bmesh.new()  # remains in ob local coordinates
    bm.from_object(ob, context.scene,
        deform=True, render=False, cage=False, face_normals=True)
    bm.faces.ensure_lookup_table()  # update bmesh index
    if matrix:
        bm.transform(matrix)
    tree = mathutils.bvhtree.BVHTree.FromBMesh(bm, epsilon=epsilon_len)
    return bm, tree

def _get_intersected_faces(bm, tree, other_tree):
    """Get intersected faces between trees."""
    overlap = tree.overlap(other_tree)
    if overlap:
        ifaces = {i_pair[0] for i_pair in overlap}
        return [bm.faces[iface] for iface in ifaces]
    return list()

def check_intersections(context, ob, other_obs=None):
    """Check ob self-intersection and intersection with other_obs."""
    DEBUG and print("BFDS: check_quality: intersections:", ob, other_obs)
    bpy.ops.object.mode_set(mode='OBJECT')
    epsilon_len = context.scene.bf_config_min_edge_length
    bad_faces = list()
    bm, tree = _get_bm_and_tree(context, ob, epsilon_len=epsilon_len)
    # Get self-intersections
    bad_faces.extend(_get_intersected_faces(bm, tree, tree))
    # Get intersections
    for other_ob in other_obs or tuple():
        matrix = ob.matrix_world.inverted() * other_ob.matrix_world
        other_bm, other_tree = _get_bm_and_tree(
            context, other_ob, epsilon_len=epsilon_len, matrix=matrix,
            )
        bad_faces.extend(_get_intersected_faces(bm, tree, other_tree))
    # Raise
    if bad_faces:
        msg = "Intersection detected, bad faces selected."
        _raise_bad_geometry(context, ob, bm, msg, bad_faces=bad_faces)          


# Raise bad geometry

def _raise_bad_geometry(context, ob, bm, msg,
        bad_verts=None, bad_edges=None, bad_faces=None):
    """Select bad elements, show them, raise BFException."""
    # Deselect all in bmesh
    for vert in bm.verts:
        vert.select = False
    bm.select_flush(False)
    # Select bad elements
    select_type = None
    if bad_faces:
        select_type = 'FACE'
        for b in bad_faces:
            b.select = True
    if bad_edges:
        select_type = 'EDGE'
        for b in bad_edges:
            b.select = True
    if bad_verts:
        select_type = 'VERT'
        for b in bad_verts:
            b.select = True
    bm.to_mesh(ob.data)
    bm.free()
    # Select object and go to edit mode
    bpy.ops.object.mode_set(mode='OBJECT')
    bpy.ops.object.select_all(action='DESELECT')
    ob.select = True
    context.scene.objects.active = ob
    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.mesh.select_mode(use_extend=False, use_expand=False, type=select_type)
    raise BFException(ob, msg)
