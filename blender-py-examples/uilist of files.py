import bpy
from bpy.types import Operator, Panel, UIList, PropertyGroup, Menu
from bpy.props import EnumProperty, CollectionProperty, IntProperty

class RENDER_UL_scene_list(UIList):
    # The draw_item function is called for each item of the collection that is visible in the list.
    #   data is the RNA object containing the collection,
    #   item is the current drawn item of the collection,
    #   icon is the "computed" icon for the item (as an integer, because some objects like materials or textures
    #   have custom icons ID, which are not available as enum items).
    #   active_data is the RNA object containing the active property for the collection (i.e. integer pointing to the
    #   active item of the collection).
    #   active_propname is the name of the active property (use 'getattr(active_data, active_propname)').
    #   index is index of the current item in the collection.
    #   flt_flag is the result of the filtering process for this item.
    #   Note: as index and flt_flag are optional arguments, you do not have to use/declare them here if you don't
    #         need them.
    def draw_item(self, context, layout, data, item, icon, active_data, active_propname):
        sce = context.scene
        ob = data
        # draw_item must handle the three layout types... Usually 'DEFAULT' and 'COMPACT' can share the same code.
        if self.layout_type in {'DEFAULT', 'COMPACT'}:
            # You should always start your row layout by a label (icon + text), or a non-embossed text field,
            # this will also make the row easily selectable in the list! The later also enables ctrl-click rename.
            # We use icon_value of label, as our given icon is an integer value, not an enum ID.
            # Note "data" names should never be translated!
                layout.prop(item, "enabled", text="")
                layout.prop_search(item, "name", bpy.data, "scenes")
        # 'GRID' layout type should be as compact as possible (typically a single icon!).
        elif self.layout_type in {'GRID'}:
            pass

# And now we can use this list everywhere in Blender. Here is a small example panel.
class RENDER_PT_scene_list(Panel):
    """Creates a Panel in the Object properties window"""
    bl_label = "Scene list"
    bl_idname = "RENDER_PT_scene_list"
    bl_space_type = 'PROPERTIES'
    bl_region_type = 'WINDOW'
    bl_context = "render"

    def draw(self, context):
        layout = self.layout

        scene = context.scene


        row = layout.row()
        # template_list now takes two new args.
        # The first one is the identifier of the registered UIList to use (if you want only the default list,
        # with no custom draw code, use "UI_UL_list").
        row.template_list("RENDER_UL_scene_list", "", scene, "render_scene_list", scene, "render_scene_list_index")
            
        col = row.column()
        sub = col.column(True)
        sub.operator(RENDER_OT_scene_list_add.bl_idname, text="", icon="ZOOMIN")
        sub.operator(RENDER_OT_scene_list_remove.bl_idname, text="", icon="ZOOMOUT")
        sub.menu(RENDER_OT_scene_list_specials.bl_idname, text="", icon="DOWNARROW_HLT")
        
        sub = col.column(True)
        sub.separator()
        sub.operator(RENDER_OT_scene_list_move.bl_idname, text="", icon="TRIA_UP").direction = 'UP'
        sub.operator(RENDER_OT_scene_list_move.bl_idname, text="", icon="TRIA_DOWN").direction = 'DOWN'


class RENDER_OT_scene_list_clear(Operator):
    bl_idname = "render.scene_list_clear"
    bl_label = "Clear"
        
    def execute(self, context):
        context.scene.render_scene_list.clear()
        return {'FINISHED'}


class RENDER_OT_scene_list_add(Operator):
    bl_idname = "render.scene_list_add"
    bl_label = "Add"
    
    def execute(self, context):
        s = context.scene
        item = s.render_scene_list.add()
        #item.name = ...
        return {'FINISHED'}


class RENDER_OT_scene_list_remove(Operator):
    bl_idname = "render.scene_list_remove"
    bl_label = "Remove"
    bl_description = "Remove item from Sprite Sheet Filelist"
    
    @classmethod
    def poll(cls, context):
        s = context.scene
        return len(s.render_scene_list) &gt; s.render_scene_list_index &gt;= 0
    
    def execute(self, context):
        s = context.scene
        s.render_scene_list.remove(s.render_scene_list_index)
        if s.render_scene_list_index &gt; 0:
            s.render_scene_list_index -= 1
        return {'FINISHED'}


class RENDER_OT_scene_list_move(Operator):
    bl_idname = "render.scene_list_move"
    bl_label = "Move"

    direction = EnumProperty(items=(
        ('UP', "Up", "Move up"),
        ('DOWN', "Down", "Move down"))
    )
    
    @classmethod
    def poll(cls, context):
        s = context.scene
        return len(s.render_scene_list) &gt; s.render_scene_list_index &gt;= 0

    def execute(self, context):
        s = context.scene
        d = -1 if self.direction == 'UP' else 1
        new_index = (s.render_scene_list_index + d) % len(s.render_scene_list)
        s.render_scene_list.move(s.render_scene_list_index, new_index)
        s.render_scene_list_index = new_index
        return {'FINISHED'}

class RENDER_OT_scene_list_specials(Menu):
    bl_idname = "RENDER_OT_scene_list_specials"
    bl_label = "Specials"

    def draw(self, context):
        layout = self.layout
        s = context.scene
        layout.label("%i Scenes" % len(s.render_scene_list), icon="SCENE_DATA")
        layout.separator()
        layout.operator(RENDER_OT_scene_list_clear.bl_idname, icon='X')

class RenderSceneEntry(bpy.types.PropertyGroup):
    # name = bpy.props.StringProperty() # this exists by default
    enabled = bpy.props.BoolProperty(name="Enabled")


def register():
    bpy.utils.register_module(__name__)
    # There's no global / singleton ID type and WindowManager doesn't serialize,
    # use Scene in lack of a better solution
    bpy.types.Scene.render_scene_list = CollectionProperty(type=RenderSceneEntry)
    bpy.types.Scene.render_scene_list_index = IntProperty(min=0)


def unregister():
    bpy.utils.register_module(__name__)
    del bpy.types.Scene.render_scene_list
    del bpy.types.Scene.render_scene_list_index


if __name__ == "__main__":
    register()



