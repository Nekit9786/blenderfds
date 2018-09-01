"""BlenderFDS, custom list operators."""

import bpy

class CUSTOM_LIST_OT():
    data_name = None   # in "scene", "object", "material", ...
    index_name = None  # Name of the active index property
    collection_name = None   # Name of the collection property

class slot_add(CUSTOM_LIST_OT):
    """Add slot to custom list."""
    bl_idname = "custom_list.slot_add"
    bl_label = "Add"
    bl_description = "Add slot"

    def set_item(self, context, item):
        """Set item in the new slot."""
        pass

    def execute(self, context):
        """Add slot and set item."""
        celem = getattr(context, self.data_name)
        cindex = getattr(celem, self.index_name)
        clist = getattr(context.scene, self.collection_name)
        item = clist.add()
        self.set_item(context, item)
        setattr(context.scene, self.index_name, len(clist)-1)
        return {"FINISHED"}

class slot_rm(CUSTOM_LIST_OT):
    """Remove slot from custom list."""
    bl_idname = "bf_catf.slot_rm"
    bl_label = "Remove"
    bl_description = "Remove slot"

    data_name = "scene"
    index_name = "bf_catf_files_index"
    collection_name = "bf_catf_files"

    @classmethod
    def poll(cls, context):
        celem = getattr(context, cls.data_name)
        return getattr(context.scene, cls.collection_name)

    def invoke(self, context, event):
        celem = getattr(context, self.data_name)
        cindex = getattr(celem, self.index_name)
        clist = getattr(context.scene, self.collection_name)
        if cindex < 0:  # available item?
            return {"FINISHED"}
        clist.remove(cindex)
        setattr(context.scene, self.index_name, max(0, cindex-1))
        return {"FINISHED"}

class slot_mv(CUSTOM_LIST_OT):
    """Move slot from custom list."""
    bl_idname = "bf_catf.slot_mv"
    bl_label = "Move"
    bl_description = "Move slot"
    direction = bpy.props.EnumProperty(items=(('UP', 'Up', ""),
                                              ('DOWN', 'Down', ""),))

    @classmethod
    def poll(cls, context):
        celem = getattr(context, cls.data_name)
        return getattr(celem, cls.collection_name)

    def execute(self, context):
        celem = getattr(context, self.data_name)
        cindex = getattr(celem, self.index_name)
        clist = getattr(context.scene, self.collection_name)
        delta = -1 if self.direction == 'UP' else 1
        neighbor = cindex + delta
        clist.move(neighbor, cindex)
        setattr(context.scene, self.index_name, max(0, min(neighbor, len(clist)-1)))
        return{'FINISHED'}

