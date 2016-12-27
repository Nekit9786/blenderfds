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

