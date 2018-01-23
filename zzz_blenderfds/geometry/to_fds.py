"""BlenderFDS, translate Blender object geometry to FDS notation."""

import bpy, bmesh
from time import time
from .geom_utils import *
from .voxelize import voxelize, pixelize
from ..exceptions import BFException

DEBUG = False

### to None

def ob_to_none(context, ob):
    DEBUG and print("BFDS: geometry.ob_to_none:", ob.name)
    return (), ""

### to XB

def ob_to_xbs_voxels(context, ob) -> "((x0,x1,y0,y1,z0,z1,), ...), 'Message'":
    """Transform ob solid geometry in XBs notation (voxelization). Never send None."""
    DEBUG and print("BFDS: geometry.ob_to_xbs_voxels:", ob.name)
    t0 = time()
    xbs, voxel_size, timing = voxelize(context, ob)
    if not xbs: return (), "No voxel created"
    msg = "{0} voxels, resolution {1:.3f} m, in {2:.0f} s".format(len(xbs), voxel_size, time()-t0)
    if DEBUG: msg += " (s:{0[0]:.0f} 1f:{0[1]:.0f}, 2g:{0[2]:.0f}, 3g:{0[3]:.0f})".format(timing)
    return xbs, msg

def ob_to_xbs_pixels(context, ob) -> "((x0,x1,y0,y1,z0,z0,), ...), 'Message'":
    """Transform ob flat geometry in XBs notation (flat voxelization). Never send None."""
    DEBUG and print("BFDS: geometry.ob_to_xbs_pixels:", ob.name)
    t0 = time()
    xbs, voxel_size, timing = pixelize(context, ob)
    if not xbs: return (), "No pixel created"
    msg = "{0} pixels, resolution {1:.3f} m, in {2:.0f} s".format(len(xbs), voxel_size, time()-t0)
    if DEBUG: msg += " (s:{0[0]:.0f} 1f:{0[1]:.0f}, 2g:{0[2]:.0f}, 3g:{0[3]:.0f})".format(timing)
    return xbs, msg

def ob_to_xbs_bbox(context, ob) -> "((x0,x1,y0,y1,z0,z1,), ...), 'Message'":
    """Transform ob solid geometry in XBs notation (bounding box). Never send None."""
    DEBUG and print("BFDS: geometry.ob_to_xbs_bbox:", ob.name)
    x0, x1, y0, y1, z0, z1 = get_global_bbox(context, ob)
    return [(x0, x1, y0, y1, z0, z1,),], ""

def ob_to_xbs_faces(context, ob) -> "((x0,x1,y0,y1,z0,z0,), ...), 'Message'": # TODO use BMesh
    """Transform ob faces in XBs notation (faces). Never send None."""
    DEBUG and print("BFDS: geometry.ob_to_xbs_faces:", ob.name)
    # Init
    result = list()
    me = get_global_mesh(context, ob)
    tessfaces = get_tessfaces(context, me)
    # For each tessface...
    for tessface in tessfaces:
        vertices = [me.vertices[vertex] for vertex in tessface.vertices]
        # Calc the bounding box in global coordinates
        bbminx, bbminy, bbminz = vertices[0].co
        bbmaxx, bbmaxy, bbmaxz = vertices[0].co
        for vertex in vertices[1:]:
            x, y, z = vertex.co
            bbminx, bbminy, bbminz = min(bbminx, x), min(bbminy, y), min(bbminz, z)
            bbmaxx, bbmaxy, bbmaxz = max(bbmaxx, x), max(bbmaxy, y), max(bbmaxz, z)
        bbd = [(bbmaxx - bbminx, 2), (bbmaxy - bbminy, 1), (bbmaxz - bbminz, 0)]
        bbd.sort()
        if bbd[0][1] == 2: bbmaxx = bbminx = (bbminx+bbmaxx)/2
        if bbd[0][1] == 1: bbmaxy = bbminy = (bbminy+bbmaxy)/2
        if bbd[0][1] == 0: bbmaxz = bbminz = (bbminz+bbmaxz)/2
        result.append((bbminx, bbmaxx, bbminy, bbmaxy, bbminz, bbmaxz,),)
    result.sort()
    # Clean up
    bpy.data.meshes.remove(me, do_unlink=True)
    # Return
    msg = len(result) > 1 and "{0} faces".format(len(result)) or ""
    return result, msg

def ob_to_xbs_edges(context, ob) -> "((x0,x1,y0,y1,z0,z1,), ...), 'Message'": # TODO use BMesh
    """Transform ob faces in XBs notation (faces). Never send None."""
    DEBUG and print("BFDS: geometry.ob_to_xbs_edges:", ob.name)
    # Init
    result = list()
    me = get_global_mesh(context, ob)
    # For each edge...
    for edge in me.edges:
        pt0x, pt0y, pt0z = me.vertices[edge.vertices[0]].co
        pt1x, pt1y, pt1z = me.vertices[edge.vertices[1]].co
        result.append((pt0x, pt1x, pt0y, pt1y, pt0z, pt1z,),)
    result.sort()
    # Clean up
    bpy.data.meshes.remove(me, do_unlink=True)
    # Return
    msg = len(result) > 1 and "{0} edges".format(len(result)) or ""
    return result, msg

# Caller function (ob.bf_xb)

choose_to_xbs = {
    "NONE"   : ob_to_none,
    "BBOX"   : ob_to_xbs_bbox,
    "VOXELS" : ob_to_xbs_voxels,
    "FACES"  : ob_to_xbs_faces,
    "PIXELS" : ob_to_xbs_pixels,
    "EDGES"  : ob_to_xbs_edges,
}

def ob_to_xbs(context, ob) -> "((x0,x1,y0,y1,z0,z0,), ...), 'Message'":
    """Transform Blender object geometry according to ob.bf_xb to FDS notation. Never send None."""
    # not ob.get("ob_to_xbs_cache") -> precalc not available or modified input conditions
    DEBUG and print("BFDS: geometry.ob_to_xbs:", ob.name)
    if not ob.get("ob_to_xbs_cache"): # ob.is_updated does not work here, checked in the handler
        ob["ob_to_xbs_cache"] = choose_to_xbs[ob.bf_xb](context, ob) # Calculate
    return ob["ob_to_xbs_cache"]

### to XYZ

def ob_to_xyzs_vertices(context, ob) -> "((x0,y0,z0,), ...), 'Message'": # TODO use BMesh
    """Transform ob vertices in XYZs notation. Never send None."""
    DEBUG and print("BFDS: geometry.ob_to_xyzs_vertices:", ob.name)
    # Init
    result = list()
    me = get_global_mesh(context, ob)
    # For each vertex...
    for vertex in me.vertices:
        pt0x, pt0y, pt0z = vertex.co
        result.append((pt0x, pt0y, pt0z,),)
    result.sort()
    # Clean up
    bpy.data.meshes.remove(me, do_unlink=True)
    # Return
    msg = len(result) > 1 and "{0} vertices".format(len(result)) or ""
    return result, msg

def ob_to_xyzs_center(context, ob) -> "((x0,y0,z0,), ...), 'Message'":
    """Transform ob center in XYZs notation. Never send None."""
    DEBUG and print("BFDS: geometry.ob_to_xyzs_center:", ob.name)
    return [(ob.location[0], ob.location[1], ob.location[2],),], ""

# Caller function (ob.bf_xyz)

choose_to_xyzs = {
    "NONE"     : ob_to_none,
    "CENTER"   : ob_to_xyzs_center,
    "VERTICES" : ob_to_xyzs_vertices,
}

def ob_to_xyzs(context, ob):
    """Transform Blender object geometry according to ob.bf_xyz to FDS notation. Never send None."""
    # not ob.get("ob_to_xyzs_cache") -> precalc not available or modified input conditions
    DEBUG and print("BFDS: geometry.ob_to_xyzs:", ob.name)
    if not ob.get("ob_to_xyzs_cache"): # ob.is_updated does not work here, checked in the handler
        ob["ob_to_xyzs_cache"] = choose_to_xyzs[ob.bf_xyz](context, ob) # Calculate
    return ob["ob_to_xyzs_cache"]
    

### to PB

def ob_to_pbs_planes(context, ob) -> "(('X',x3,), ('X',x7,), ('Y',y9,), ...), 'Message'":
    """Transform ob faces in PBs notation. Never send None."""
    DEBUG and print("BFDS: geometry.ob_to_pbs_planes:", ob.name)
    # Init
    result = list()
    xbs, msg = ob_to_xbs_faces(context, ob)
    # For each face build a plane...
    for xb in xbs:
        if   abs(xb[1] - xb[0]) < epsilon: result.append((0,xb[0],),) # PBX is 0
        elif abs(xb[3] - xb[2]) < epsilon: result.append((1,xb[2],),) # PBY is 1
        elif abs(xb[5] - xb[4]) < epsilon: result.append((2,xb[4],),) # PBZ is 2
        else: raise ValueError("BFDS: Building planes impossible, problem in ob_to_xbs_faces.")
    result.sort()
    # Nothing to clean up, return
    msg = len(result) > 1 and "{0} planes".format(len(result)) or ""
    return result, msg

# Caller function (ob.bf_pb)

choose_to_pbs = {
    "NONE"   : ob_to_none,
    "PLANES" : ob_to_pbs_planes,
}

def ob_to_pbs(context, ob):
    """Transform Blender object geometry according to ob.bf_pb to FDS notation. Never send None."""
    # not ob.get("ob_to_pbs_cache") -> precalc not available or modified input conditions
    DEBUG and print("BFDS: geometry.ob_to_pbs:", ob.name)
    if not ob.get("ob_to_pbs_cache"): # ob.is_updated does not work here, checked in the handler
        ob["ob_to_pbs_cache"] = choose_to_pbs[ob.bf_pb](context, ob) # Calculate
    return ob["ob_to_pbs_cache"]

### to GEOM

# Current format:
# &GEOM ID='FEM_MESH',
#       SURF_ID='CONE',
#       MATL_ID='CONE',
#       VERTS=0.0699,-0.0146,3.4286,
#             0.0714, 0.0000,3.4286,
#             ...
#       FACES=1,2,3,
#             4,5,2,
#             ...   
#       /

# Future format:
# &GEOM ID='FEM_MESH',
#       SURF_IDV='CONE','Gypsum plaster','Wood',
#       MATL_ID='Matl1',
#       VERTS=0.0699,-0.0146,3.4286,
#             0.0714, 0.0000,3.4286,
#             ...
#       FACES=1,2,3, 1,
#             4,5,2, 1,
#             ...
#       /
#
# Notes:
# ob.data.materials[me.polygons[2].material_index]
# ob.data.materials[me.polygons[0].material_index]
# ob.materials?
#

def ob_to_geom(context, ob) -> "verts, faces":
    """Transform Blender object geometry to GEOM FDS notation. Never send a None."""
    assert(ob.type == 'MESH')
    epsilon = .000001 # FIXME global epsilon
    
    # Get the new bmesh from the Object, apply modifiers, set in global coordinates, and triangulate
    bm = bmesh.new()
    bm.from_object(ob, context.scene, deform=True, render=False, cage=False, face_normals=True)
    bm.transform(ob.matrix_world)
    bmesh.ops.triangulate(bm, faces=bm.faces)

    # Check self intersection
    import mathutils
    tree = mathutils.bvhtree.BVHTree.FromBMesh(bm, epsilon=0.00001) # FIXME epsilon
    if tree.overlap(tree): raise BFException(ob, "Object self intersection detected.")

    # Check edges:
    # - manifold, each edge should join two faces, no more no less
    # - contiguous normals, adjoining faces should have normals in the same directions
    # - no degenerate edges, zero lenght edges
    for edge in bm.edges:
        if not edge.is_manifold: raise BFException(ob, "Non manifold edges detected.")
        if not edge.is_contiguous: raise BFException(ob, "Adjoining faces have opposite normals.")
        if edge.calc_length() <= epsilon: raise BFException(ob, "Zero lenght edges detected.")
    
    # Check degenerate faces, zero area faces
    for face in bm.faces:
        if face.calc_area() <= epsilon: raise BFException(ob, "Zero area faces detected.")

    # Check loose vertices, vertices that have no connectivity
    for vert in bm.verts:
        if not bool(vert.link_edges): raise BFException(ob, "Loose vertices detected.")

    # Get its vertex coordinates
    verts = [(v.co.x, v.co.y, v.co.z) for v in bm.verts]

    # Get its faces by vertex index
    faces = [(f.verts[0].index+1, f.verts[1].index+1, f.verts[2].index+1) for f in bm.faces]

    # Free the bmesh
    bm.free()

    # Set up msg
    msg = "{} vertices, {} faces".format(len(verts), len(faces))
            
    return verts, faces, msg
