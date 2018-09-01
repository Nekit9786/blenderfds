"""BlenderFDS, algorithms for triangulated surfaces."""

import bpy, bmesh, mathutils
from time import time
from math import floor, ceil

from ..exceptions import BFException
from . import utils

DEBUG = True


def get_trisurface(context, ob) -> "mas, verts, faces":
    """Get triangulated surface from object ready for FDS GEOM format."""
    # Check and init
    DEBUG and print("BFDS: get_triangles")
    assert(ob.type == 'MESH')
    if not ob.data.vertices:
        raise BFException(ob, "Empty object!")
    # Check original mesh quality
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
        raise BFException(ob, "No referenced SURF, add at least one Material.")
    for material_slot in material_slots:
        ma = material_slot.material
        if not ma.bf_export:
            bpy.data.objects.remove(ob_tmp, True)
            raise BFException(ob, "Referenced SURF ID='{}' is not exported.".format(ma.name))
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
    msg = "Non manifold or open geometry detected, bad edges selected."
    _raise_bad_geometry(context, ob, bm, msg, bad_edges=bad_edges)
    # Check manifold vertices
    for vert in bm.verts:
        if not vert.is_manifold:
            bad_verts.append(vert)
    msg = "Non manifold vertices detected, bad vertices selected."
    _raise_bad_geometry(context, ob, bm, msg, bad_verts=bad_verts)
    # Check contiguous normals, adjoining faces should have normals
    # in the same directions
    for edge in bm.edges:
        if not edge.is_contiguous:
            bad_edges.append(edge)
    msg = "Inconsistent face normals detected, bad edges selected."
    _raise_bad_geometry(context, ob, bm, msg, bad_edges=bad_edges)
    # Check no degenerate edges, zero lenght edges
    for edge in bm.edges:
        if edge.calc_length() <= epsilon_len:
            bad_edges.append(edge)
    msg = "Too short edges detected, bad edges selected."
    _raise_bad_geometry(context, ob, bm, msg, bad_edges=bad_edges)
    # Check degenerate faces, zero area faces
    for face in bm.faces:
        if face.calc_area() <= epsilon_area:
            bad_faces.append(face)
    msg = "Too small area faces detected, bad faces selected."
    _raise_bad_geometry(context, ob, bm, msg, bad_faces=bad_faces)
    # Check loose vertices, vertices that have no connectivity
    for vert in bm.verts:
        if not bool(vert.link_edges):
            bad_verts.append(vert)
    msg = "Loose vertices detected, bad vertices selected."
    _raise_bad_geometry(context, ob, bm, msg, bad_verts=bad_verts)
    # Check duplicate vertices
    _check_duplicate_vertices(context, ob, bm, epsilon_len)
    # Check inverted normals
    original_normal = bm.faces[0].normal.copy()
    bmesh.ops.recalc_face_normals(bm, faces=bm.faces)  # recalc outside
    if original_normal.dot(bm.faces[0].normal) < 0.:  # opposed normals?
        raise BFException(ob, "Face normals are pointing towards the inside, "
            "update needed.")
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
    msg = "Duplicate vertices detected, bad vertices selected."
    _raise_bad_geometry(context, ob, bm, msg, bad_verts=bad_verts)

def check_intersections(obs, context):  # FIXME test and make operator
    """Check self and mutual intersection of objects."""
    DEBUG and print("BFDS: check_quality: intersection:", obs)
    for i, ob in enumerate(obs):
        # Check self intersection
        epsilon_len = context.scene.bf_config_min_edge_length
        tree = mathutils.bvhtree.BVHTree.FromObject(context, ob.scene, epsilon=epsilon_len)
        _check_intersection(context, ob, tree, epsilon_len, other_ob=ob, other_tree=tree)
        # Check intersection with others
        for other_ob in obs[i+1:]:
            _check_intersection(context, ob, tree, epsilon_len, other_ob, other_tree=None)

def _check_intersection(context, ob, tree, epsilon_len, other_ob, other_tree=None):
    """Check intersection between ob and other_ob."""
    DEBUG and print("BFDS: check_quality: intersection:", ob.name, other_ob.name)
    bad_faces = list()
    # First check bboxes for intersection (fast)
    bb0 = [ob.matrix_world * mathutils.Vector(b) for b in ob.bound_box]
    bb1 = [other_ob.matrix_world * mathutils.Vector(b) for b in other_ob.bound_box]
    if bb0[6][0] < bb1[0][0] or bb1[6][0] < bb0[0][0] or \
       bb0[6][1] < bb1[0][1] or bb1[6][1] < bb0[0][1] or \
       bb0[6][2] < bb1[0][2] or bb1[6][2] < bb0[0][2]:
        DEBUG and print("BFDS: check_quality: intersection: none (fast): ", ob.name, other_ob.name)
        return
    # Then check tree (slow)
    if not other_tree:
        other_tree = mathutils.bvhtree.BVHTree.FromObject(
            other_ob, deform=True, render=True, cage=True, epsilon=epsilon_len)
    overlap = tree.overlap(other_tree)
    if overlap:
        ifaces = {i_pair[0] for i_pair in overlap}  # only from first tree
        bad_faces = [bm.faces[iface] for iface in ifaces]
        msg = "Intersection detected, bad faces selected."
        _raise_bad_geometry(context, ob, bm, msg, bad_faces=bad_faces)
    DEBUG and print("BFDS: check_quality: intersection: none (slow): ", ob.name, other_ob.name)

def _raise_bad_geometry(context, ob, bm, msg,
        bad_verts=None, bad_edges=None, bad_faces=None):
    """Select bad elements, show them, raise BFException."""
    if not bad_verts and not bad_edges and not bad_faces:
        return
    # Deselect faces, edges, verts in bmesh
    for face in bm.faces:
        face.select = False
    for edge in bm.edges:
        edge.select = False
    for vert in bm.verts:
        vert.select = False
    # Select bad elements
    if bad_faces:
        for b in bad_faces:
            b.select = True
    if bad_edges:
        for b in bad_edges:
            b.select = True
    if bad_verts:
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
    bpy.ops.mesh.select_mode(use_extend=False, use_expand=False, type='VERT')
    raise BFException(ob, msg)
