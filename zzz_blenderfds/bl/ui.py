"""BlenderFDS, menus and other ui mods."""

import bpy
from bpy.types import Panel, Header
from bpy.utils import register_class

from . import operators

DEBUG = False


# Register/Unregister

def register():
    """Register menus and other ui mods."""
    # Properties for UI simplification
    bpy.types.WindowManager.bf_sp_context = bpy.props.EnumProperty(
        items=_sp_items,
        update=_sp_items_update,
        default="OBJECT"
    )
    # Simplify UI, if requested from user's preferences
    if bpy.context.user_preferences.addons["zzz_blenderfds"].preferences.bf_pref_simplify_ui:
        # Simplify info editor (upper menu)
        from .simplified_ui import space_info
        for cls in space_info.classes:
            try: register_class(cls)
            except ValueError: pass
        # Simplify view3d editor
        from .simplified_ui import space_view3d
        for cls in space_view3d.classes:
            try: register_class(cls)
            except ValueError: pass
        # Simplify space properties header
        from .simplified_ui import space_properties
        for cls in space_properties.classes:
            try: register_class(cls)
            except ValueError: pass
        # Simplify modifiers
        from .simplified_ui import properties_data_modifier
        for cls in properties_data_modifier.classes:
            try: register_class(cls)
            except ValueError: pass
        # Treat (rewire or unregister) unused Blender bpy.types
        _treat_unused_bl_classes()

    # Append import/export menus
    bpy.types.INFO_MT_file_export.prepend(operators.export_OT_fds_case_menu)
    bpy.types.INFO_MT_file_import.prepend(operators.import_OT_fds_snippet_menu)
    bpy.types.INFO_MT_file_import.prepend(operators.import_OT_fds_case_menu)

def unregister():
    """Unregister menus and other ui mods"""
    # info menu
    bpy.types.INFO_MT_editor_menus.draw_menus = _INFO_MT_editor_menus_draw_menus_tmp


# Rewire draw functions

def _INFO_MT_editor_menus_draw_menus_tmp(layout, context):
    """Show restart needed on top"""
    row = layout.row()
    row.alert = True
    row.operator("wm.quit_blender", text="Restart neeeded", icon='QUIT')


def _VIEW3D_PT_tools_add_object_draw(self, context):  # FIXME insert into recoded UI
    layout = self.layout
    col = layout.column(align=True)
    self.draw_add_mesh(col, label=True)


def _unused_header_draw(self, context):
    """Draw generic unused header."""
    layout = self.layout
    row = layout.row(align=True)
    row.template_header()


# Rewire space properties header

_sp_items = (
    ('SCENE','Scene','Scene','SCENE_DATA',1),
    ('OBJECT','Object','Object','OBJECT_DATA',3),
    ('MATERIAL','Material','Material','MATERIAL',5),
    ('MODIFIER','Modifiers','Object modifiers','MODIFIER',10),
)


def _sp_items_update(self, context):
    # Get the right space on screen (Properties Panel)
    space = context.space_data
    # When not called by the Properties Panel
    if space and space.type != 'PROPERTIES':
        print("_sp_items_update Not called by Properties Panel")  # FIXME
        space = None
        for window in context.window_manager.windows:
            for area in window.screen.areas:
                if area.type == 'PROPERTIES':
                    space = area.spaces[0]
                    break
    if not space:
        return
    # Update Properties Panel context
    try:
        space.context = self.bf_sp_context
    except TypeError:
        self.bf_sp_context = 'SCENE'


# Treat (rewire or unregister) unused Blender bpy.types

# Configuration

_used_bl_space_type = (
    "TEXTEDITOR", "TEXT_EDITOR", "USER_PREFERENCES", "CONSOLE",
    "FILE_BROWSER", "INFO", "OUTLINER",
)

_used_headers = (
    'Header', 'PROPERTIES_HT_header', 'VIEW3D_HT_header',
)

_used_panels = (
    'Panel',
    'SCENE_PT_BF_HEAD', 'SCENE_PT_BF_TIME', 'SCENE_PT_BF_DUMP', 'SCENE_PT_BF_MISC', 'SCENE_PT_BF_REAC',
    'OBJECT_PT_context_object',
    'OBJECT_PT_BF_MESH', 'OBJECT_PT_BF_EMPTY', 'OBJECT_PT_BF_TMP',
    'MATERIAL_PT_BF',
    'MATERIAL_PT_context_material', 'MaterialButtonsPanel',
    'Cycles_PT_context_material', 'CyclesButtonsPanel',
    'DATA_PT_modifiers', 'RENDER_PT_render',
    'OBJECT_PT_relations',
)

_unused_panels = (
        "VIEW3D_PT_view3d_shading", "VIEW3D_PT_view3d_motion_tracking",
        "VIEW3D_PT_view3d_name",  # "VIEW3D_PT_transform_orientations",
        "VIEW3D_PT_context_properties",
        "VIEW3D_PT_tools_meshweight",
        "VIEW3D_PT_grease_pencil",
    )

_used_panel_by_bl_space_type = "PROPERTIES", "VIEW_3D"
# Next was: "Create", "Relations", "Options", "Grease Pencil"
_used_panel_by_bl_category = "Tools",
_used_panel_by_bl_region_type = "UI"


# Treat unused classes

def _treat_unused_bl_classes():
    """Treat (rewire or unregister) unused Blender bpy.types ."""
    for bt_name in dir(bpy.types):
        # Init
        bt_cls = getattr(bpy.types, bt_name, None)
        bt_bl_space_type = getattr(bt_cls, "bl_space_type", None)
        # Surely used bts
        if bt_bl_space_type in _used_bl_space_type:
            continue
        # Other Headers and Panels
        # Panels
        if issubclass(bt_cls, Panel):
            if bt_name in _used_panels:
                continue
            if bt_name not in _unused_panels and bt_bl_space_type in _used_panel_by_bl_space_type:
                bt_bl_category = getattr(bt_cls, "bl_category", None)
                if bt_bl_category and bt_bl_category in _used_panel_by_bl_category:
                    continue
                bt_bl_region_type = getattr(bt_cls, "bl_region_type", None)
                if bt_bl_region_type and bt_bl_region_type in _used_panel_by_bl_region_type:
                    continue
            # If nothing else applies, unregister the Panel
            if DEBUG:
                print("BFDS: Unregister Panel:", bt_name)
            bpy.utils.unregister_class(bt_cls)
            continue
        # Headers
        if issubclass(bt_cls, Header):
            if bt_name in _used_headers:
                continue
            if DEBUG:
                print("BFDS: Rewire Header:", bt_name)
            # Unused header, rewire its draw function
            bt_cls.draw = _unused_header_draw
            continue
