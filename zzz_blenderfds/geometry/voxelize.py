"""BlenderFDS, voxelize algorithm."""

import bpy
from time import time

import bmesh # FIXME
from math import floor, ceil # FIXME

from ..types import BFException
from .geom_utils import * 
from . import tmp_objects

DEBUG = False

# TODO port to numpy for speed

# FIXME not used, not finished
# FIXME I wish I could merge flat pixels and faces!

def pixelize(context, ob) -> "(xbs, voxel_size, timing)":
    """Pixelize object."""
    print("BFDS: voxelize.pixelize:", ob.name)

    ## Init: check, voxel_size
    t0 = time()
    if not ob.data.vertices: raise BFException(ob, "Empty object!")
    if ob.bf_xb_custom_voxel: voxel_size = ob.bf_xb_voxel_size
    else: voxel_size = context.scene.bf_default_voxel_size

    # Get a copy of original object in global coordinates
    me_bvox = get_global_mesh(context, ob)
    ob_bvox = get_new_object(context, context.scene, "bvox", me_bvox, linked=False)

    # Get bbox in global coordinates
    bbox = get_bbox(ob_bvox)

    # Launch infinite rays parallel to x axis in object bbox
    # Trasform intercepted faces in square faces perp to x axis
    # Then grow them in z and y directions, finally transform to xbs
    # Same for y
    # Same for z

# Geometry Utilities (mathutils.geometry) or raycast
#mathutils.geometry.intersect_ray_tri(v1, v2, v3, ray, orig, clip=True)
#Returns the intersection between a ray and a triangle, if possible, returns None otherwise.

#Parameters:	
#v1 (mathutils.Vector) – Point1
#v2 (mathutils.Vector) – Point2
#v3 (mathutils.Vector) – Point3
#ray (mathutils.Vector) – Direction of the projection
#orig (mathutils.Vector) – Origin
#clip (boolean) – When False, don’t restrict the intersection to the area of the triangle, use the infinite plane defined by the triangle.
#Returns:	
#The point of intersection or None if no intersection is found
#Return type:	
#mathutils.Vector or None


# "global" coordinates are absolute coordinate referring to Blender main origin of axes,
# that are directly transformed to FDS coordinates (that refer to the only origin of axes) 

def voxelize(context, ob, flat=False) -> "(xbs, voxel_size, timing)":
    """Voxelize object."""
    print("BFDS: voxelize.voxelize:", ob.name)
    
    # ob: original object in local coordinates
    # ob_cvox: original object in global coordinates with normalized bounding box
    # ob_bvox: original object in global coordinates, before voxelization
    # ob_avox: voxelized object in global coordinates, after voxelization

    ## Init: check, get voxel_size
    t0 = time()
    if not ob.data.vertices: raise BFException(ob, "Empty object!")
    voxel_size = _get_voxel_size(ob)

    ## Normalize object global bbox FIXME mi piacerebbe usare sempre lo stesso oggetto, non doverne creare tre!
    # The best bbox: has an odd number of voxels,
    # is at least 1 voxel far from the object,
    # is aligned with the origin, as FDS MESHes
    me_cvox = get_global_mesh(context, ob)
    ob_cvox = get_new_object(context, context.scene, "cvox", me_cvox, linked=False)
    # Get bbox
    bbox = get_bbox(ob_cvox)
    # Alignment and minimum distance
    pv0 = ( # in voxels
        floor(bbox[0]/voxel_size) - 1,
        floor(bbox[2]/voxel_size) - 1,
        floor(bbox[4]/voxel_size) - 1,
    )
    pv1 = ( # in voxels
        ceil(bbox[1]/voxel_size) + 1,
        ceil(bbox[3]/voxel_size) + 1,
        ceil(bbox[5]/voxel_size) + 1,
    )
    # Correct pv1 for odd number of voxels
    pv1 = ( # in voxels
        pv0[0] + (pv1[0]-pv0[0])//2*2+1,
        pv0[1] + (pv1[1]-pv0[1])//2*2+1,
        pv0[2] + (pv1[2]-pv0[2])//2*2+1,
    )
    # New bbox and vertices
    bbox = ( # in meters
        pv0[0] * voxel_size, pv1[0] * voxel_size,
        pv0[1] * voxel_size, pv1[1] * voxel_size,
        pv0[2] * voxel_size, pv1[2] * voxel_size,
    )
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
    # Insert vertices into mesh # FIXME spostare in geom_utils.py
    bm = bmesh.new()
    # convert the current mesh to a bmesh (must be in edit mode)
    #bpy.ops.object.mode_set(mode='EDIT')
    bm.from_mesh(me_cvox)
    #bpy.ops.object.mode_set(mode='OBJECT')  # return to object mode
    for v in verts: bm.verts.new(v)  # add a new vert
    # make the bmesh the object's mesh
    bm.to_mesh(me_cvox)  
    bm.free()  # always do this when finished

    ## Voxelize object
    # Get a copy of original object in global coordinates (remesh works in local coordinates)
    me_bvox = get_global_mesh(context, ob_cvox) # FIXME ob_cvox instead of ob
    ob_bvox = get_new_object(context, context.scene, "bvox", me_bvox, linked=False)
    # If flat, solidify and get flatten function for later generated xbs
    if flat: flat_origin, choose_flatten = _solidify_flat_ob(context, ob_bvox, voxel_size/3.)
    # Apply remesh modifier to ob_bvox
    octree_depth, scale = _calc_remesh_modifier(context, ob_bvox, voxel_size) # ob_bvox not ob!
    _apply_remesh_modifier(context, ob_bvox, octree_depth, scale)
    # Get voxelized object
    ob_avox = get_new_object(context, context.scene, "avox", get_global_mesh(context, ob_bvox), linked=False)

    ## Find, build and grow boxes
    # Get and check tessfaces
    tessfaces = get_tessfaces(context, ob_avox.data)
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
    # If flat, flatten xbs at flat_origin
    if flat: xbs = choose_flatten(xbs, flat_origin)

    ## Clean up
    if DEBUG:
#        bpy.data.objects.remove(ob_cvox) # FIXME levare        
#        bpy.data.objects.remove(ob_bvox)        
        context.scene.objects.link(ob_cvox) # create unlinked for speed # FIXME levare
        tmp_objects.tmp_set(context, ob, ob_cvox) # it is left as a tmp object # FIXME levare
        context.scene.objects.link(ob_bvox) # create unlinked for speed # FIXME levare
        tmp_objects.tmp_set(context, ob, ob_bvox) # it is left as a tmp object # FIXME levare

        context.scene.objects.link(ob_avox) # create unlinked for speed
        tmp_objects.tmp_set(context, ob, ob_avox) # it is left as a tmp object
    else:
        bpy.data.objects.remove(ob_bvox)
        bpy.data.objects.remove(ob_avox)

    ## Return
    return xbs, voxel_size, (t2-t1, t4-t3, t5-t4, t6-t5) # this is timing: sort, 1b, 2g, 3g 

def _get_voxel_size(ob) -> "voxel_size":
    """Get object voxel size."""
    if ob.bf_xb_custom_voxel: return ob.bf_xb_voxel_size
    else: return context.scene.bf_default_voxel_size

def _solidify_flat_ob(context, ob, thickness):
    """Solidify a flat object. Apply modifier, return flat_origin and flatten function for later generated xbs."""
    # Set flat_origin at any vertices of original flat ob
    flat_origin = ob.data.vertices[0].co
    # Choose flat dimension and set according xbs flatten function
    if   ob.dimensions[0] < epsilon: choose_flatten = _x_flatten_xbs # the face is normal to x axis
    elif ob.dimensions[1] < epsilon: choose_flatten = _y_flatten_xbs # ... to y axis
    elif ob.dimensions[2] < epsilon: choose_flatten = _z_flatten_xbs # ... to z axis
    else: raise BFException(ob, "Not flat and normal to axis, cannot create pixels.")
    # Solidify
    _apply_solidify_modifier(context, ob, thickness)
    # Return the function that is going to be used to flatten later generated xbs
    return flat_origin, choose_flatten

# When appling a remesh modifier, object max dimension is scaled up FIXME explain better
# and divided in 2 ** octree_depth voxels
# Example: dimension = 4.2, voxel_size = 0.2, octree_depth = 5, number of voxels = 2^5 = 32, scale = 3/4 = 0.75
# |-----|-----|-----|-----| voxels, dimension / scale
#    |=====.=====.=====|    dimension

def _calc_remesh_modifier(context, ob_bvox, voxel_size):
    """Calc Remesh modifier parameters for voxel_size."""
    # Get max dimension and init flag
    dimension = max(ob_bvox.dimensions)
    dimension_too_large = True
    # Find righ octree_depth and relative scale
    for octree_depth in range(1,13):
        scale = dimension / voxel_size / 2 ** octree_depth
        if 0.010 < scale < 0.990:
            dimension_too_large = False
            break
    if dimension_too_large: raise BFException(ob_bvox, "Too large for desired resolution, split object!")
    if DEBUG:
        print("Remesh: dimension, voxel_size, 2 ** octree_depth, scale:", dimension, voxel_size, 2 ** octree_depth, scale) # FIXME
    # Return
    return octree_depth, scale

def _apply_remesh_modifier(context, ob, octree_depth, scale):
    """Apply remesh modifier for voxelization."""
    mo = ob.modifiers.new('voxels_tmp','REMESH') # apply modifier
    mo.mode, mo.use_remove_disconnected, mo.octree_depth, mo.scale = 'BLOCKS', False, octree_depth, scale

def _apply_solidify_modifier(context, ob, thickness):
    """Apply solidify modifier with centered thickness."""
    mo = ob.modifiers.new('solid_tmp','SOLIDIFY') # apply modifier
    mo.thickness, mo.offset = thickness, 0. # Set centered thickness # FIXME check

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

def _x_flatten_xbs(xbs, flat_origin) -> "[(l0, l0, y0, y1, z0, z1), ...]":
    """Flatten voxels to obtain pixels (normal to x axis) at flat_origin height."""
    print("BFDS: _x_flatten_xbs:", len(xbs))
    return [[flat_origin[0], flat_origin[0], xb[2], xb[3], xb[4], xb[5]] for xb in xbs]
        
def _y_flatten_xbs(xbs, flat_origin) -> "[(x0, x1, l0, l0, z0, z1), ...]":
    """Flatten voxels to obtain pixels (normal to x axis) at flat_origin height."""
    print("BFDS: _y_flatten_xbs:", len(xbs))
    return [[xb[0], xb[1], flat_origin[1], flat_origin[1], xb[4], xb[5]] for xb in xbs]

def _z_flatten_xbs(xbs, flat_origin) -> "[(x0, x1, y0, y1, l0, l0), ...]":
    """Flatten voxels to obtain pixels (normal to x axis) at flat_origin height."""
    print("BFDS: _z_flatten_xbs:", len(xbs))
    return [[xb[0], xb[1], xb[2], xb[3], flat_origin[2], flat_origin[2]] for xb in xbs]
