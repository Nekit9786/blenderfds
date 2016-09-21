"""BlenderFDS, temporary object management."""

import bpy

def restore_all(context):
    """Restore all original obs, delete all tmp objects, delete all cached geometry."""
    if context.mode != 'OBJECT': bpy.ops.object.mode_set(mode='OBJECT', toggle=False)
    sc = context.scene
    # For all objects in this scene
    for ob in sc.objects:
        # Unlink and delete temporary object
        if ob.bf_is_tmp:
            sc.objects.unlink(ob)
            bpy.data.objects.remove(ob)
            continue
        # Restore original object and delete cached geometry
        if ob.bf_has_tmp: ob.bf_has_tmp, ob.hide = False, False
        ob["ob_to_xbs_cache"] = False
        ob["ob_to_xyzs_cache"] = False
        ob["ob_to_pbs_cache"] = False

def del_my_tmp(context, ob):
    """Restore original ob, delete my tmp object."""
    if context.mode != 'OBJECT': bpy.ops.object.mode_set(mode='OBJECT', toggle=False)
    sc = context.scene
    # Restore original object
    if ob.bf_has_tmp: ob.bf_has_tmp, ob.hide = False, False
    # Unlink and delete my temporary object
    for child in ob.children:
        if child.bf_is_tmp:
            sc.objects.unlink(child)
            bpy.data.objects.remove(child)
            break
    # Set ob as selected and active object
    ob.select = True
    sc.objects.active = ob

def tmp_set(context, ob, ob_tmp):
    """Link ob_tmp as temporary object of ob."""
    # Set original object
    ob.bf_has_tmp = True
    ob.hide = True
    # Set temporary object
    ob_tmp.bf_is_tmp = True
    ob_tmp.active_material = ob.active_material
    ob_tmp.layers = ob.layers
    ob_tmp.show_wire = True
    ob_tmp.select = True
    # Set parenting and keep position
    ob_tmp.parent = ob
    ob_tmp.matrix_parent_inverse = ob.matrix_world.inverted()
