"""BlenderFDS, voxelize algorithm."""

import bpy
from time import time

from math import floor, ceil

from ..types import BFException
from .geom_utils import * 
from . import tmp_objects

DEBUG = False

# TODO port to numpy for speed
# TODO port to bmesh

# "global" coordinates are absolute coordinate referring to Blender main origin of axes,
# that are directly transformed to FDS coordinates (that refers its coordinates to the
# one and only origin of axes) 

def pixelize(context, ob) -> "(xbs, voxel_size, timing)":
    """Pixelize object."""
    print("BFDS: voxelize.pixelize:", ob.name)
    # Init
    flat_axis = _get_flat_axis(ob)
    voxel_size = _get_voxel_size(context, ob)
    thickness = voxel_size
    bbox = get_global_bbox(context,ob)
    flat_origin = bbox[0],bbox[2],bbox[4]
    # Solidify
    choose_flatten = (_x_flatten_xbs, _y_flatten_xbs, _z_flatten_xbs)[flat_axis]
    ob_solid = _get_solidify_ob(context, ob, thickness)
    ob_solid.bf_xb_voxel_size = voxel_size # prepare new object, you are not passing the voxel size
    ob_solid.bf_xb_custom_voxel = True
    # Voxelize
    xbs, voxel_size, ts = voxelize(context, ob_solid)
    # Flatten
    xbs = choose_flatten(xbs, flat_origin)
    ## Clean up
    if DEBUG: ob_solid.set_tmp(context, ob)
    else: bpy.data.objects.remove(ob_solid, do_unlink=True)
    # Return
    return xbs, voxel_size, ts

def _get_solidify_ob(context, ob, thickness) -> "ob":
    """Get a new unlinked obj with solidify modifier for voxelization applied."""
    ob_new = get_new_object(
        context,
        context.scene,
        "solidified_tmp",
        me=ob.data,
        linked=False
    )
    # Create modifier
    mo = ob_new.modifiers.new('solidify_tmp','SOLIDIFY')
    mo.thickness = thickness
    mo.offset = 0. # centered
    # Apply modifier
    me = ob_new.to_mesh(scene=context.scene, apply_modifiers=True, settings="RENDER")
    ob_new.data = me
    # Return
    return ob_new

def voxelize(context, ob) -> "(xbs, voxel_size, timing)":
    """Voxelize object."""
    print("BFDS: voxelize.voxelize:", ob.name)
    # Init
    t0 = time()
    voxel_size = _get_voxel_size(context, ob)
    # Check
    if not ob.data.vertices: raise BFException(ob, "Empty object!")
    ## Prepare ob
    ob_norm = _get_normalized_ob(context, ob, voxel_size)
    ob_remesh = _get_remesh_ob(context, ob_norm, voxel_size)
    ## Find, build and grow boxes
    # Get and check tessfaces
    tessfaces = get_tessfaces(context, ob_remesh.data)
    if not tessfaces: raise BFException(ob, "No tessfaces available, cannot voxelize.")
    # Sort tessfaces centers by face normal: normal to x, to y, to z.
    t1 = time()
    x_tessfaces, y_tessfaces, z_tessfaces = _sort_tessfaces_by_normal(tessfaces)
    # Choose fastest procedure: less tessfaces => less time required
    # Better use the smallest collection first!
    t2 = time()
    choose = [
        (len(x_tessfaces), x_tessfaces, _x_tessfaces_to_boxes, _grow_boxes_along_x, _x_boxes_to_xbs),
        (len(y_tessfaces), y_tessfaces, _y_tessfaces_to_boxes, _grow_boxes_along_y, _y_boxes_to_xbs),
        (len(z_tessfaces), z_tessfaces, _z_tessfaces_to_boxes, _grow_boxes_along_z, _z_boxes_to_xbs),
    ]
    choose.sort(key=lambda k:k[0]) # sort by len(tessfaces)
    # Build minimal boxes along 1st axis, using floors
    t3 = time()
    boxes, origin = choose[0][2](choose[0][1], voxel_size) # eg. _x_tessfaces_to_boxes(x_tessfaces, voxel_size)
    # Grow boxes along 2nd axis
    t4 = time()
    boxes = choose[1][3](boxes) # eg. _grow_boxes_along_y(boxes)
    # Grow boxes along 3rd axis
    t5 = time()
    boxes = choose[2][3](boxes) # eg. _grow_boxes_along_z(boxes)
    ## Make xbs
    # Transform grown boxes in xbs
    t6 = time()
    xbs = choose[0][4](boxes, voxel_size, origin) # eg. _x_boxes_to_xbs(boxes, ...)
    ## Clean up
    if DEBUG:
        ob_norm.set_tmp(context, ob)
        ob_remesh.set_tmp(context, ob)
    else:
        bpy.data.objects.remove(ob_norm, do_unlink=True)
        bpy.data.objects.remove(ob_remesh, do_unlink=True)
    ## Return
    return xbs, voxel_size, (t2-t1, t4-t3, t5-t4, t6-t5) # this is timing: sort, 1b, 2g, 3g 

def _get_voxel_size(context, ob) -> "voxel_size":
    """Get voxel size for object."""
    if ob.bf_xb_custom_voxel: return ob.bf_xb_voxel_size
    else: return context.scene.bf_default_voxel_size

def _get_normalized_ob(context, ob, voxel_size) -> "ob":
    """Get a new, global, unlinked object with normalized bbox."""
    # Get global mesh (absolute origin, all modifiers applied) and create a new object    
    me_norm = get_global_mesh(context, ob)
    ob_norm = get_new_object(
        context,
        context.scene,
        "normalized_tmp",
        me_norm,
        linked=False
    )
    # Get its bbox (it's already global) and an origin
    bbox = get_bbox(ob_norm)
    # Calc vertices for voxel alignment to origin (as FDS does)
    # and minimum distance (in voxel unit)
    pv0 = (
        floor(bbox[0]/voxel_size)-1, # at least one voxel far from obj
        floor(bbox[2]/voxel_size)-1,
        floor(bbox[4]/voxel_size)-1,
    )
    pv1 = (
        ceil(bbox[1]/voxel_size)+1,
        ceil(bbox[3]/voxel_size)+1,
        ceil(bbox[5]/voxel_size)+1,
    )
    # Correct pv1 for odd number of voxels dimension (in voxel unit)
    pv1 = (
        pv0[0] + (pv1[0]-pv0[0])//2*2+1,  # at least one voxel far from obj
        pv0[1] + (pv1[1]-pv0[1])//2*2+1,
        pv0[2] + (pv1[2]-pv0[2])//2*2+1,
    )
    # Calc new bbox (in Blender units)
    bbox = (
        pv0[0] * voxel_size, pv1[0] * voxel_size,
        pv0[1] * voxel_size, pv1[1] * voxel_size,
        pv0[2] * voxel_size, pv1[2] * voxel_size,
    )
    # Setup and insert new vertices
    verts = (
        (bbox[0], bbox[2], bbox[4]),
        (bbox[0], bbox[2], bbox[5]),
        (bbox[0], bbox[3], bbox[4]),
        (bbox[0], bbox[3], bbox[5]),
        (bbox[1], bbox[2], bbox[4]),
        (bbox[1], bbox[2], bbox[5]),
        (bbox[1], bbox[3], bbox[4]),
        (bbox[1], bbox[3], bbox[5]),
    )
    insert_vertices_into_mesh(me_norm, verts)
    # Update mesh and return
    ob_norm.data = me_norm
    return ob_norm

# When appling a remesh box modifier, object max dimension is scaled up FIXME explain better VERY BAD!
# and subdivided in 2 ** octree_depth voxels - 1 cubic voxels 
# Example: dimension = 4.2, voxel_size = 0.2, octree_depth = 5, number of voxels = 2^5-1 = 31, scale = 3/4 = 0.75
# |-----|-----|-----|-----| voxels, dimension / scale
#    |=====.=====.=====|    dimension

def _get_remesh_ob(context, ob, voxel_size) -> "ob":
    """Get a new unlinked obj with remesh modifier for voxelization applied."""
    # Init modifier
    dimension = max(ob.dimensions)
    octree_depth, scale = _calc_remesh(voxel_size, dimension)
    # Create modifier
    mo = ob.modifiers.new('remesh_tmp','REMESH')
    mo.mode = 'BLOCKS'
    mo.use_remove_disconnected = False
    mo.octree_depth = octree_depth
    mo.scale = scale
    # Apply modifier
    me = ob.to_mesh(scene=context.scene, apply_modifiers=True, settings="RENDER")
    return bpy.data.objects.new("remeshed_tmp", me)

def _calc_remesh(voxel_size, dimension):
    """Calc Remesh modifier parameters for voxel_size and dimension."""
    octree_depth = 0.
    while True:
        octree_depth += 1.
        scale = dimension / voxel_size / 2 ** octree_depth
        if 0.010 < scale < 0.990: break
    if DEBUG:
        print("BFDS: _calc_remesh_modifier:")
        print("    dimension:", dimension)
        print("    voxel_size:", voxel_size)
        print("    octree_depth:", octree_depth)
        print("    scale:", scale)
    return octree_depth, scale

# Sort tessfaces by normal: collection of tessfaces normal to x, to y, to z
# tessfaces created by the Remesh modifier in BLOCKS mode are perpendicular to a local axis
# we used a global object, the trick is done: tessfaces are perpendicular to global axis

def _sort_tessfaces_by_normal(tessfaces):
    """Sort tessfaces: normal to x axis, y axis, z axis."""
    print("BFDS: _sort_tessfaces_by_normal:", len(tessfaces))
    x_tessfaces, y_tessfaces, z_tessfaces = list(), list(), list()
    for tessface in tessfaces:
        normal = tessface.normal
        if   abs(normal[0]) > .9: x_tessfaces.append(tessface) # tessface is normal to x axis
        elif abs(normal[1]) > .9: y_tessfaces.append(tessface) # ... to y axis
        elif abs(normal[2]) > .9: z_tessfaces.append(tessface) # ... to z axis
        else: raise ValueError("BFDS: voxelize._sort_tessfaces_by_normal: abnormal face")
    return x_tessfaces, y_tessfaces, z_tessfaces

# First, we transform the global tessface center coordinates
# in integer coordinates referred to origin point:
# - voxel_size is used as step;
# - origin is the first of tessface centers;
# - (center[0] - origin[0]) / voxel_size is rounded from float to integers: ix, iy, iz

# Then we pile integer heights of floors for each location:
# (ix, iy -> location int coordinates):
#    (iz0, iz1, ... -> list of floors int coordinates)

# Last we use this "floor levels" (eg. izs) to detect solid volumes.
# Eg. at location (ix, iy) of int coordinates, at izs[0] floor go into solid,
# at izs[1] go out of solid, at izs[2] go into solid, ...
# z axis --> floor 0|==solid==1| void 2|==solid==3| void ...
# If solid is manifold, len(izs) is an even number: go into solid at izs[0], get at last out of it at izs[-1].

# In fact this floors can be easily transformed in boxes:
# (ix0, ix1, iy0, iy1, iz0, iz1)
# boxes are very alike XBs, but in integer coordinates.

def _x_tessfaces_to_boxes(x_tessfaces, voxel_size) -> "[(ix0, ix1, iy0, iy1, iz0, iz1), ...], origin":
    """Transform _x_tessfaces into minimal boxes."""
    print("BFDS: _x_tessfaces_to_boxes:", len(x_tessfaces))
    # Create floors
    origin = tuple(x_tessfaces[0].center) # First tessface center is origin
    floors = dict() # {(3,4):(3,4,15,25,), (3,5):(3,4,15,25), ...}
    for tessface in x_tessfaces:
        center = tuple(tessface.center)
        ix = round((center[0] - origin[0]) / voxel_size) # integer coordinates for this face (a floor)
        iy = round((center[1] - origin[1]) / voxel_size)
        iz = round((center[2] - origin[2]) / voxel_size)
        try: floors[(iy, iz)].append(ix) # append face ix to list of ixs
        except KeyError: floors[(iy, iz)] = [ix,] # or create new list of ixs from face ix
    # Create minimal boxes
    boxes = list()
    while floors:
        (iy, iz), ixs = floors.popitem()
        ixs.sort() # sort from bottom to top in +x direction
        while ixs:
            ix1 = ixs.pop() # pop from top to bottom in -x direction
            ix0 = ixs.pop()
            boxes.append((ix0, ix1, iy, iy, iz, iz,))
    return boxes, origin

def _y_tessfaces_to_boxes(y_tessfaces, voxel_size) -> "[(ix0, ix1, iy0, iy1, iz0, iz1), ...], origin":
    """Transform _y_tessfaces into minimal boxes."""
    print("BFDS: _y_tessfaces_to_boxes:", len(y_tessfaces))
    # Create floors
    origin = tuple(y_tessfaces[0].center) # First tessface center becomes origin
    floors = dict() # {(3,4):(3,4,15,25,), (3,5):(3,4,15,25), ...}
    for tessface in y_tessfaces:
        center = tuple(tessface.center)
        ix = round((center[0] - origin[0]) / voxel_size) # integer coordinates for this face (a floor)
        iy = round((center[1] - origin[1]) / voxel_size)
        iz = round((center[2] - origin[2]) / voxel_size)
        try: floors[(ix, iz)].append(iy) # append face iy to existing list of iys
        except KeyError: floors[(ix, iz)] = [iy,] # or create new list of iys from face iy
    # Create minimal boxes
    boxes = list()
    while floors:
        (ix, iz), iys = floors.popitem()
        iys.sort() # sort from bottom to top in +y direction
        while iys:
            iy1 = iys.pop() # pop from top to bottom in -y direction
            iy0 = iys.pop()
            boxes.append((ix, ix, iy0, iy1, iz, iz,))
    return boxes, origin

def _z_tessfaces_to_boxes(z_tessfaces, voxel_size) -> "[(ix0, ix1, iy0, iy1, iz0, iz1), ...], origin":
    """Transform _z_tessfaces into minimal boxes."""
    print("BFDS: _z_tessfaces_to_boxes:", len(z_tessfaces))
    # Create floors
    origin = tuple(z_tessfaces[0].center) # First tessface center becomes origin
    floors = dict() # {(3,4):(3,4,15,25,), (3,5):(3,4,15,25), ...}
    for tessface in z_tessfaces:
        center = tuple(tessface.center)
        ix = round((center[0] - origin[0]) / voxel_size) # integer coordinates for this face (a floor)
        iy = round((center[1] - origin[1]) / voxel_size)
        iz = round((center[2] - origin[2]) / voxel_size)
        try: floors[(ix, iy)].append(iz) # append face iz to existing list of izs
        except KeyError: floors[(ix, iy)] = [iz,] # or create new list of izs from face iz
    # Create minimal boxes  
    boxes = list()
    while floors:
        (ix, iy), izs = floors.popitem()
        izs.sort() # sort from bottom to top in +z direction
        while izs:
            iz1 = izs.pop() # pop from top to bottom in -z direction
            iz0 = izs.pop()
            boxes.append((ix, ix, iy, iy, iz0, iz1,))
    return boxes, origin

# Merge each minimal box with available neighbour boxes in axis direction

def _grow_boxes_along_x(boxes) -> "[(ix0, ix1, iy0, iy1, iz0, iz1), ...]":
    """Grow boxes by merging neighbours along x axis."""
    print("BFDS: _grow_boxes_along_x:", len(boxes))
    boxes_grown = list()
    while boxes:
        ix0, ix1, iy0, iy1, iz0, iz1 = boxes.pop()
        while True: # grow into +x direction
            box_desired = (ix1+1, ix1+1, iy0, iy1, iz0, iz1,)
            try: boxes.remove(box_desired)
            except ValueError: break
            ix1 += 1
        while True: # grow into -x direction
            box_desired = (ix0 - 1, ix0 - 1, iy0, iy1, iz0, iz1,)
            try: boxes.remove(box_desired)
            except ValueError: break
            ix0 -= 1
        boxes_grown.append((ix0, ix1, iy0, iy1, iz0, iz1))
    return boxes_grown

def _grow_boxes_along_y(boxes) -> "[(ix0, ix1, iy0, iy1, iz0, iz1), ...]":
    """Grow boxes by merging neighbours along y axis."""
    print("BFDS: _grow_boxes_along_y:", len(boxes))
    boxes_grown = list()
    while boxes:
        ix0, ix1, iy0, iy1, iz0, iz1 = boxes.pop()
        while True: # grow into +y direction
            box_desired = (ix0, ix1, iy1+1, iy1+1, iz0, iz1)
            try: boxes.remove(box_desired)
            except ValueError: break
            iy1 += 1
        while True: # grow into -y direction
            box_desired = (ix0, ix1, iy0 - 1, iy0 - 1, iz0, iz1)
            try: boxes.remove(box_desired)
            except ValueError: break
            iy0 -= 1
        boxes_grown.append((ix0, ix1, iy0, iy1, iz0, iz1))
    return boxes_grown

def _grow_boxes_along_z(boxes) -> "[(ix0, ix1, iy0, iy1, iz0, iz1), ...]":
    """Grow boxes by merging neighbours along z axis."""
    print("BFDS: _grow_boxes_along_z:", len(boxes))
    boxes_grown = list()
    while boxes:
        ix0, ix1, iy0, iy1, iz0, iz1 = boxes.pop()
        while True: # grow into +z direction
            box_desired = (ix0, ix1, iy0, iy1, iz1+1, iz1+1)
            try: boxes.remove(box_desired)
            except ValueError: break
            iz1 += 1
        while True: # grow into -z direction
            box_desired = (ix0, ix1, iy0, iy1, iz0 - 1, iz0 - 1)
            try: boxes.remove(box_desired)
            except ValueError: break
            iz0 -= 1
        boxes_grown.append((ix0, ix1, iy0, iy1, iz0, iz1))
    return boxes_grown

# Trasform boxes in int coordinates to xbs in global coordinates

def _x_boxes_to_xbs(boxes, voxel_size, origin) -> "[(x0, x1, y0, y1, z0, z1), ...]":
    """Trasform boxes (int coordinates) to xbs (global coordinates)."""
    # Init
    print("BFDS: _x_boxes_to_xbs:", len(boxes))
    xbs = list()
    voxel_size_half = voxel_size / 2. + epsilon # This epsilon is used to obtain overlapping boxes
    # Build xbs
    while boxes:
        ix0, ix1, iy0, iy1, iz0, iz1 = boxes.pop()
        x0, y0, z0 = ( # origin + location + movement to lower left corner
            origin[0] + ix0 * voxel_size - epsilon, # - 0., this is already at floor level
            origin[1] + iy0 * voxel_size - voxel_size_half,
            origin[2] + iz0 * voxel_size - voxel_size_half,
        )
        x1, y1, z1 = ( # origin + location + movement to upper right corner
            origin[0] + ix1 * voxel_size + epsilon, # + 0., this is already at floor level
            origin[1] + iy1 * voxel_size + voxel_size_half,
            origin[2] + iz1 * voxel_size + voxel_size_half,
        )
        xbs.append([x0, x1, y0, y1, z0, z1],)
    return xbs

def _y_boxes_to_xbs(boxes, voxel_size, origin) -> "[(x0, x1, y0, y1, z0, z1), ...]":
    """Trasform boxes (int coordinates) to xbs (global coordinates)."""
    print("BFDS: _y_boxes_to_xbs:", len(boxes))
    # Init
    xbs = list()
    voxel_size_half = voxel_size / 2. + epsilon # This epsilon is used to obtain overlapping boxes
    # Build xbs
    while boxes:
        ix0, ix1, iy0, iy1, iz0, iz1 = boxes.pop()
        x0, y0, z0 = ( # origin + location + movement to lower left corner
            origin[0] + ix0 * voxel_size - voxel_size_half,
            origin[1] + iy0 * voxel_size - epsilon, # - 0., this is already at floor level
            origin[2] + iz0 * voxel_size - voxel_size_half,
        )
        x1, y1, z1 = ( # origin + location + movement to upper right corner
            origin[0] + ix1 * voxel_size + voxel_size_half,
            origin[1] + iy1 * voxel_size + epsilon, # + 0., this is already at floor level
            origin[2] + iz1 * voxel_size + voxel_size_half,
        )
        xbs.append([x0, x1, y0, y1, z0, z1],)
    return xbs

def _z_boxes_to_xbs(boxes, voxel_size, origin) -> "[(x0, x1, y0, y1, z0, z1), ...]":
    """Trasform boxes (int coordinates) to xbs (global coordinates)."""
    # Init
    print("BFDS: _z_boxes_to_xbs:", len(boxes))
    xbs = list()
    voxel_size_half = voxel_size / 2. + epsilon # This epsilon is used to obtain overlapping boxes
    # Build xbs
    while boxes:
        ix0, ix1, iy0, iy1, iz0, iz1 = boxes.pop()
        x0, y0, z0 = ( # origin + location + movement to lower left corner
            origin[0] + ix0 * voxel_size - voxel_size_half,
            origin[1] + iy0 * voxel_size - voxel_size_half,
            origin[2] + iz0 * voxel_size - epsilon, # - 0., this is already at floor level
        )
        x1, y1, z1 = ( # origin + location + movement to upper right corner
            origin[0] + ix1 * voxel_size + voxel_size_half,
            origin[1] + iy1 * voxel_size + voxel_size_half,
            origin[2] + iz1 * voxel_size + epsilon, # + 0., this is already at floor level
        )
        xbs.append([x0, x1, y0, y1, z0, z1],)
    return xbs

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

