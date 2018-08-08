"""BlenderFDS, voxelization algorithms."""

import bpy, bmesh
from time import time
from math import floor, ceil

from ..exceptions import BFException
from . import utils
from . import tmp_objects

DEBUG = True


def get_voxels(context, ob, align_voxels=True):
    """Get voxels from object in xbs format."""
    # Check and init
    DEBUG and print("BFDS: calc_voxels.get_voxels")
    t0 = time()
    assert(ob.type == 'MESH')
    if not ob.data.vertices:
        raise BFException(ob, "Empty object!")
    voxel_size = _get_voxel_size(context, ob)
    # Create new object, and link it
    ob_tmp = utils.object_get_global_copy(context, ob, suffix='_vox_tmp')
    # Align voxels to global origin
    if align_voxels:
        _align_remesh_to_global_origin(context, ob_tmp, voxel_size)
    # Create remesh modifier
    octree_depth, scale = _init_remesh_mod(context, ob_tmp, voxel_size)
    mo = ob_tmp.modifiers.new('remesh_tmp','REMESH')
    mo.mode, mo.use_remove_disconnected = 'BLOCKS', False
    mo.octree_depth, mo.scale = octree_depth, scale
    # Apply remesh modifier and remove it
    ob_tmp.data = ob_tmp.to_mesh(
        scene=context.scene,
        apply_modifiers=True,
        settings="RENDER",
        calc_tessface=True,
        calc_undeformed=False,
    )
    ob_tmp.modifiers.remove(mo)
    # Get faces and sort them according to normals
    t1 = time()
    faces = ob_tmp.data.tessfaces
    x_faces, y_faces, z_faces = list(), list(), list()
    for face in faces:
        normal = face.normal
        if   abs(normal[0]) > .9:
            x_faces.append(face)  # face is normal to x axis
        elif abs(normal[1]) > .9:
            y_faces.append(face)  # ... to y axis
        elif abs(normal[2]) > .9:
            z_faces.append(face)  # ... to z axis
        else:
            raise ValueError("BFDS: abnormal face")
    # Choose shorter list of faces, relative functions, and parameters
    t2 = time()
    choices = [
        (len(x_faces), _get_boxes_along_x, x_faces, _grow_boxes_along_x, 0),
        (len(y_faces), _get_boxes_along_y, y_faces, _grow_boxes_along_y, 2),
        (len(z_faces), _get_boxes_along_z, z_faces, _grow_boxes_along_z, 4),
    ]
    choices.sort(key=lambda choice: choice[0])
    get_boxes = choices[0][1]  # get boxes according to fastest orientation
    faces = choices[0][2]
    grow_boxes_along_first_axis = choices[1][3]   # 1st axis growing direction
    first_sort_by = choices[2][4]
    grow_boxes_along_second_axis = choices[2][3]  # 2nd axis growing direction
    second_sort_by = choices[1][4]
    # For each face find other sides and build boxes data structure
    t3 = time()
    boxes, origin = get_boxes(faces, voxel_size)
    # Join boxes along other axis and return their global coordinates
    t4 = time()
    boxes = grow_boxes_along_first_axis(boxes, first_sort_by)
    t5 = time()
    boxes = grow_boxes_along_second_axis(boxes, second_sort_by)
    # Clean up
    t6 = time()
    bpy.data.objects.remove(ob_tmp, do_unlink=True)
    xbs = list(_get_box_xbs(boxes, origin, voxel_size))
    # Return with timing: sort, 1b, 2g, 3g
    return xbs, voxel_size, (t2-t1, t4-t3, t5-t4, t6-t5)

# When appling a remesh modifier to a Blender Object in BLOCK mode,
# the object max dimension is scaled up and divided in
# 2 ** octree_depth voxels - 1 cubic voxels
# Example: dimension = 4.2, voxel_size = 0.2,
# octree_depth = 5, number of voxels = 2^5-1 = 31,
# scale = 3/4 = 0.75
# The next function reverses the procedures and calculate octree_depth
# and scale that generate the desired voxel_size.

def _init_remesh_mod(context, ob, voxel_size) -> "octree_depth, scale":
    """Calc remesh modifier."""
    dimension = max(ob.dimensions)
    octree_depth = 0.
    while True:
        octree_depth += 1.
        scale = dimension / voxel_size / 2 ** octree_depth
        if 0.010 < scale < 0.990:
            break
        if octree_depth > 9:
            raise BFException(ob, "Object too large for its voxel size, split in parts.")
    return octree_depth, scale

def _get_voxel_size(context, ob) -> "voxel_size":
    """Get voxel_size for object."""
    if ob.bf_xb_custom_voxel:
        return ob.bf_xb_voxel_size
    else:
        return context.scene.bf_default_voxel_size

# When appling a remesh modifier to a Blender Object, the octree is aligned with
# the max dimension of the object. By inserting some loose vertices to the
# temporary object, we can align the voxelization to FDS global origin

def _align_remesh_to_global_origin(context, ob, voxel_size):
    """Modify object mesh for remesh voxel alignment to global origin."""
    bb = utils.get_bbox(ob)  # the object is already global
    # Calc new bbox (in Blender units)
    #      +---+ pv1
    #      |   |
    #  pv0 +---+
    pv0 = (
        floor(bb[0]/voxel_size)-1, # at least one voxel away from obj
        floor(bb[2]/voxel_size)-1,
        floor(bb[4]/voxel_size)-1,
    )
    pv1 = (
        ceil(bb[1]/voxel_size)+1,
        ceil(bb[3]/voxel_size)+1,
        ceil(bb[5]/voxel_size)+1,
    )
    pv1 = (
        pv0[0] + (pv1[0]-pv0[0])//2*2+1,  # at least one voxel away from obj
        pv0[1] + (pv1[1]-pv0[1])//2*2+1,  # and odd number of voxels
        pv0[2] + (pv1[2]-pv0[2])//2*2+1,
    )
    bb = (
        pv0[0] * voxel_size, pv1[0] * voxel_size,
        pv0[1] * voxel_size, pv1[1] * voxel_size,
        pv0[2] * voxel_size, pv1[2] * voxel_size,
    )
    verts = (
        (bb[0], bb[2], bb[4]), (bb[0], bb[2], bb[5]),
        (bb[0], bb[3], bb[4]), (bb[0], bb[3], bb[5]),
        (bb[1], bb[2], bb[4]), (bb[1], bb[2], bb[5]),
        (bb[1], bb[3], bb[4]), (bb[1], bb[3], bb[5]),
    )
    # Insert new vertices into ob mesh
    bm = bmesh.new()
    bm.from_mesh(ob.data)
    for v in verts:
        bm.verts.new(v)
    bm.to_mesh(ob.data)
    bm.free()
    # Push an update (circumvent bug in ob.dimensions)
    context.scene.update()

# The following functions transform remesh modifier faces into boxes,
# by raytracing along the requested axis. Each face is transformed into
# integer coordinates according to a local origin (the first face center).
# The faces are piled up, sorted, and transformed into solids:
# Eg.: z axis --> pile 0|==solid==1| void 2|==solid==3| void ...
# If solid is manifold, len(izs) is an even number:
# go into solid at izs[0], get at last out of it at izs[-1].
# In fact this piles can be easily transformed in boxes:
# (ix0, ix1, iy0, iy1, iz0, iz1)
# boxes are very alike XBs, but in integer coordinates.

# For example, with y faces:
#  y ^
#    |     B
#  3 +---+=x=+---+ F = (.5, 0): face origin (integer coordinates)
#    |   H   H   | face A = (1, 0, 0)
#  2 +---+---+---+ face B = (1, 3, 0)
#    |   H   H   | O = (0, 0): box origin (integer coordinate)
#  1 +---+---+---+ box AB = (1, 2, 0, 3, 0, 0)
#    |   H   H   |
#  0 +-F-+=x=+---+->
#    0   1 A 2   3 x

def _get_boxes_along_x(faces, voxel_size) -> "boxes, origin":
    """Get minimal boxes from faces by raytracing along x axis."""
    DEBUG and print("BFDS: _get_boxes_along_x")
    # First face center becomes origin of the integer grid for faces
    f_origin = tuple(faces[0].center)
    hvs = voxel_size / 2.
    origin = (f_origin[0], f_origin[1]-hvs, f_origin[2]-hvs)
    # Get integer coordinates of faces and
    # classify faces in integer piles along z
    # piles = {(3,4):(3,4,15,25,), (3,5):(3,4,15,25), ...}
    piles = dict()
    for face in faces:
        center = face.center
        ix, iy, iz = (
            round((center[0]-f_origin[0]) / voxel_size),
            round((center[1]-f_origin[1]) / voxel_size),
            round((center[2]-f_origin[2]) / voxel_size),
        )
        try:
            piles[(iy, iz)].append(ix)
        except KeyError:
            piles[(iy, iz)] = [ix,]
    # Create boxes by raytracing piles along axis
    # boxes = [[ix0, ix1, iy0, iy1, iz0, iz1], ...]
    boxes = list()
    for (iy, iz), ixs in piles.items():
        ixs.sort()  # sort in +x direction
        while ixs:
            ix1, ix0 = ixs.pop(), ixs.pop()  # pop solid volumes from top to bottom
            boxes.append([ix0, ix1, iy, iy+1, iz, iz+1,])
    return boxes, origin

def _get_boxes_along_y(faces, voxel_size) -> "boxes, origin":
    """Get minimal boxes from faces by raytracing along y axis."""
    DEBUG and print("BFDS: _get_boxes_along_y")
    # First face center becomes origin of the integer grid for faces
    f_origin = tuple(faces[0].center)
    hvs = voxel_size / 2.
    origin = (f_origin[0]-hvs, f_origin[1], f_origin[2]-hvs)
    # Get integer coordinates of faces and
    # classify faces in integer piles along z
    # piles = {(3,4):(3,4,15,25,), (3,5):(3,4,15,25), ...}
    piles = dict()
    for face in faces:
        center = face.center
        ix, iy, iz = (
            round((center[0]-f_origin[0]) / voxel_size),
            round((center[1]-f_origin[1]) / voxel_size),
            round((center[2]-f_origin[2]) / voxel_size),
        )
        try:
            piles[(iz, ix)].append(iy)
        except KeyError:
            piles[(iz, ix)] = [iy,]
    # Create boxes by raytracing piles along axis
    # boxes = [[ix0, ix1, iy0, iy1, iz0, iz1], ...]
    boxes = list()
    for (iz, ix), iys in piles.items():
        iys.sort()  # sort in +y direction
        while iys:
            iy1, iy0 = iys.pop(), iys.pop()  # pop solid volumes from top to bottom
            boxes.append([ix, ix+1, iy0, iy1, iz, iz+1,])
    return boxes, origin

def _get_boxes_along_z(faces, voxel_size) -> "boxes, origin":
    """Get minimal boxes from faces by raytracing along z axis."""
    DEBUG and print("BFDS: _get_boxes_along_z")
    # First face center becomes origin of the integer grid for faces
    f_origin = tuple(faces[0].center)
    hvs = voxel_size / 2.
    origin = (f_origin[0]-hvs, f_origin[1]-hvs, f_origin[2])
    # Get integer coordinates of faces and
    # classify faces in integer piles along z
    # piles = {(3,4):(3,4,15,25,), (3,5):(3,4,15,25), ...}
    piles = dict()
    for face in faces:
        center = face.center
        ix, iy, iz = (
            round((center[0]-f_origin[0]) / voxel_size),
            round((center[1]-f_origin[1]) / voxel_size),
            round((center[2]-f_origin[2]) / voxel_size),
        )
        try:
            piles[(ix, iy)].append(iz)
        except KeyError:
            piles[(ix, iy)] = [iz,]
    # Create boxes by raytracing piles along axis
    # boxes = [[ix0, ix1, iy0, iy1, iz0, iz1], ...]
    boxes = list()
    for (ix, iy), izs in piles.items():
        izs.sort()  # sort in +z direction
        while izs:
            iz1, iz0 = izs.pop(), izs.pop()  # pop solid volumes from top to bottom
            boxes.append([ix, ix+1, iy, iy+1, iz0, iz1,])
    return boxes, origin

# The following functions reduce the number of boxes in xbs format,
# used to describe the geometry, by merging them

def _grow_boxes_along_x(boxes, sort_by):
    """Grow boxes by merging neighbours along x axis."""
    DEBUG and print("BFDS: _grow_boxes_along_x")
    # Sort boxes
    boxes.sort(key=lambda box: (box[sort_by], box[0]))
    # Grow boxes in -x direction, starting from last one
    boxes_grown = list()
    box = boxes.pop()
    while boxes:
        abox = boxes.pop()
        # Check same iz0, iz1, iy0, iy1, and touching abox ix1 with box ix0
        if abox[4] == box[4] and abox[5] == box[5] and \
           abox[2] == box[2] and abox[3] == box[3] and \
           abox[1] == box[0]:
            box[0] = abox[0]  # grow box along -x
        else:
            boxes_grown.append(box)  # stash the resulting box
            box = abox  # init next cycle
    # Stash the last one
    boxes_grown.append(box)
    return boxes_grown

def _grow_boxes_along_y(boxes, sort_by):
    """Grow boxes by merging neighbours along y axis."""
    DEBUG and print("BFDS: _grow_boxes_along_y")
    # Sort boxes
    boxes.sort(key=lambda box: (box[sort_by], box[2]))
    # Grow boxes in -y direction, starting from last one
    boxes_grown = list()
    box = boxes.pop()
    while boxes:
        abox = boxes.pop()
        # Check same iz0, iz1, ix0, ix1, and touching abox iy1 with box iy0
        if abox[4] == box[4] and abox[5] == box[5] and \
           abox[0] == box[0] and abox[1] == box[1] and \
           abox[3] == box[2]:
            box[2] = abox[2]  # grow box along -y
        else:
            boxes_grown.append(box)  # stash the resulting box
            box = abox  # init next cycle
    # Stash the last one
    boxes_grown.append(box)
    return boxes_grown

def _grow_boxes_along_z(boxes, sort_by):
    """Grow boxes by merging neighbours along z axis."""
    DEBUG and print("BFDS: _grow_boxes_along_z")
    # Sort boxes
    boxes.sort(key=lambda box: (box[sort_by], box[4]))
    # Grow boxes in -z direction, starting from last one
    boxes_grown = list()
    box = boxes.pop()
    while boxes:
        abox = boxes.pop()
        # Check same iy0, iy1, ix0, ix1, and touching abox iz1 with box iz0
        if abox[2] == box[2] and abox[3] == box[3] and \
           abox[0] == box[0] and abox[1] == box[1] and \
           abox[5] == box[4]:
            box[4] = abox[4]  # grow box along -z
        else:
            boxes_grown.append(box)  # stash the resulting box
            box = abox  # init next cycle
    # Stash the last one
    boxes_grown.append(box)
    return boxes_grown

# Transform boxes in integer coordinates, back to global coordinates

def _get_box_xbs(boxes, origin, voxel_size) -> "xbs":
    """Transform boxes to xbs in global coordinates."""
    epsilon = 1E-5
    return ((
            origin[0] + box[0] * voxel_size - epsilon,
            origin[0] + box[1] * voxel_size + epsilon,
            origin[1] + box[2] * voxel_size - epsilon,
            origin[1] + box[3] * voxel_size + epsilon,
            origin[2] + box[4] * voxel_size - epsilon,
            origin[2] + box[5] * voxel_size + epsilon,
            ) for box in boxes)

# FIXME Verification case: box with 3 different sides in any orientation.
# Should always come out an only box


# "global" coordinates are absolute coordinate referring to Blender main origin of axes,
# that are directly transformed to FDS coordinates (that refers its coordinates to the
# one and only origin of axes)

def pixelize(context, ob) -> "(xbs, voxel_size, timing)":
    """Pixelize object."""
    print("BFDS: voxelize.pixelize:", ob.name)
    # Init
    voxel_size = _get_voxel_size(context, ob)
    ob_global = get_new_object(
        context,
        context.scene,
        "global_ob",
        get_global_mesh(context, ob),
        linked=False,
    )
    flat_axis = _get_flat_axis(ob_global)
    bbox = get_global_bbox(context,ob_global)
    flat_origin = bbox[0],bbox[2],bbox[4]
    choose_flatten = (_x_flatten_xbs, _y_flatten_xbs, _z_flatten_xbs)[flat_axis]
    # Solidify
    ob_solid = _get_solidify_ob(context, ob_global, voxel_size * 1.5)
    # Voxelize
    ob_solid.bf_xb_voxel_size = voxel_size # prepare new object, you are not passing the voxel size
    ob_solid.bf_xb_custom_voxel = True
    xbs, voxel_size, ts = voxelize(context, ob_solid)
    # Flatten
    xbs = choose_flatten(xbs, flat_origin)
    ## Clean up
    if DEBUG:
        ob_global.set_tmp(context, ob)
        ob_solid.set_tmp(context, ob)
    else:
        bpy.data.objects.remove(ob_global, do_unlink=True)
        bpy.data.objects.remove(ob_solid, do_unlink=True)
    # Return
    return xbs, voxel_size, ts

def _get_solidify_ob(context, ob, thickness) -> "ob":
    """Get a new unlinked obj with solidify modifier for voxelization applied."""
    ob_new = get_new_object(
        context,
        context.scene,
        "solidified_tmp",
        me=ob.data,
        linked=False,
    )
    # Create modifier
    mo = ob_new.modifiers.new('solidify_tmp','SOLIDIFY')
    mo.thickness = thickness
    mo.offset = 0. # centered
    ob_new.data.update(calc_edges=True, calc_tessface=True) # Or it will not update the mesh
    # Apply modifier
    me = ob_new.to_mesh(scene=context.scene, apply_modifiers=True, settings="RENDER")
    ob_new.data = me
    # Return
    return ob_new




# Flatten xbs to obtain pixels

def _get_flat_axis(ob):
    """Get object flat axis."""
    dimensions = ob.dimensions
    choose = [
        (dimensions[0],0), # object faces are normal to x axis
        (dimensions[1],1), # ... to y axis
        (dimensions[2],2), # ... to z axis
    ]
    choose.sort(key=lambda k:k[0]) # sort by dimension
    return choose[0][1]

def _x_flatten_xbs(xbs, flat_origin) -> "[(l0, l0, y0, y1, z0, z1), ...]":
    """Flatten voxels to obtain pixels (normal to x axis) at flat_origin height."""
    print("BFDS: _x_flatten_xbs:", len(xbs))
    return [[flat_origin[0], flat_origin[0], xb[2], xb[3], xb[4], xb[5]] for xb in xbs]

def _y_flatten_xbs(xbs, flat_origin) -> "[(x0, x1, l0, l0, z0, z1), ...]":
    """Flatten voxels to obtain pixels (normal to y axis) at flat_origin height."""
    print("BFDS: _y_flatten_xbs:", len(xbs))
    return [[xb[0], xb[1], flat_origin[1], flat_origin[1], xb[4], xb[5]] for xb in xbs]

def _z_flatten_xbs(xbs, flat_origin) -> "[(x0, x1, y0, y1, l0, l0), ...]":
    """Flatten voxels to obtain pixels (normal to z axis) at flat_origin height."""
    print("BFDS: _z_flatten_xbs:", len(xbs))
    return [[xb[0], xb[1], xb[2], xb[3], flat_origin[2], flat_origin[2]] for xb in xbs]

