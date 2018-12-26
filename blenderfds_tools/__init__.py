#    BlenderFDS, an open tool for the NIST Fire Dynamics Simulator
#    Copyright (C) 2013  Emanuele Gissi, http://www.blenderfds.org
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""BlenderFDS Tools"""

import subprocess, os, bpy
from bpy.types import Operator, AddonPreferences, Panel
from bpy.props import StringProperty, IntProperty, BoolProperty, EnumProperty, FloatProperty

bl_info = {
    "name": "BlenderFDS Tools",
    "author": "Emanuele Gissi",
    "version": (0, 0, 1),
    "blender": (2, 7, 9),
    "api": 35622,
    "location": "View3D > FDS Tools",
    "description": "Tools for BlenderFDS",
    "warning": "",
    "wiki_url": "http://www.blenderfds.org/",
    "tracker_url": "https://github.com/firetools/blenderfds/issues",
    "support": "COMMUNITY",
    "category": "3D View",
}

class BFTPref(bpy.types.AddonPreferences):
    bl_idname = __name__

    bft_im_filepath = StringProperty(
        name="InstantMeshes executable filepath",
        subtype='FILE_PATH',
        )

    def draw(self, context):
            self.layout.prop(self, "bft_im_filepath")
            self.layout.label("Download from https://github.com/wjakob/instant-meshes")


class BFTClean(bpy.types.Operator):  # FIXME
    bl_idname="ops.bft_clean"
    bl_label="Clean"
    bl_description="Clean up of selected Objects"
    bl_options={'REGISTER', 'UNDO'}
    bl_region_type="WINDOW"

    def execute(self, context):
        # Check and init
        if not bpy.context.selected_objects:
            self.report({'ERROR'}, "No selected Objects")
            return {'CANCELLED'}
        wm = context.window_manager
        bpy.ops.object.mode_set(mode='EDIT')
        bpy.ops.mesh.select_all(action='SELECT')
        bpy.ops.mesh.delete_loose()
        bpy.ops.mesh.dissolve_degenerate(threshold = wm.bft_edge)
        bpy.ops.mesh.remove_doubles(threshold = wm.bft_edge)
        return {'FINISHED'}


class BFTBlShrinkwrap(bpy.types.Operator):  # FIXME
    bl_idname="ops.bft_bl_shrinkwrap"
    bl_label="Shrinkwrap"
    bl_description="Retopology of selected Objects with Blender Shrinkwrap modifier from the bbox"
    bl_options={'REGISTER', 'UNDO'}
    bl_region_type="WINDOW"
    
    operation = bpy.props.StringProperty(name="operation", default='autorun')

    def execute(self, context):
        # Check and init
        if not bpy.context.selected_objects:
            self.report({'ERROR'}, "No selected Objects")
            return {'CANCELLED'}
        # FIXME
        # Copy original properties to imported Object
        bpy.ops.object.bf_props_to_sel_obs()
        return {'FINISHED'}
        
class BFTBlRemesh(bpy.types.Operator):
    bl_idname="ops.bft_bl_remesh"
    bl_label="Remesh"
    bl_description="Retopology of selected Objects with Blender Remesh modifier"
    bl_options={'REGISTER', 'UNDO'}
    bl_region_type="WINDOW"
    
    operation = bpy.props.StringProperty(name="operation", default='autorun')

    def execute(self, context):
        # Check and init
        if not bpy.context.selected_objects:
            self.report({'ERROR'}, "No selected Objects")
            return {'CANCELLED'}
        # Add modifier
        wm = context.window_manager
        for ob in bpy.context.selected_objects:
            octree_depth, scale = self._init_remesh_mod(context, ob, wm.bft_edge)
            mo = ob.modifiers.new('bl_remesh','REMESH')
            mo.mode, mo.use_remove_disconnected = 'SHARP', True
            mo.octree_depth, mo.scale = octree_depth, scale
        return {'FINISHED'}

    # When appling a remesh modifier to a Blender Object in BLOCK mode,
    # the object max dimension is scaled up and divided in
    # 2 ** octree_depth voxels - 1 cubic voxels
    # Example: dimension = 4.2, voxel_size = 0.2,
    # octree_depth = 5, number of voxels = 2^5-1 = 31,
    # scale = 3/4 = 0.75
    # The next function reverses the procedures and calculate octree_depth
    # and scale that generate the desired voxel_size.

    def _init_remesh_mod(self, context, ob, edge_size) -> "octree_depth, scale":
        """Calc remesh modifier."""
        dimension = max(ob.dimensions)
        octree_depth = 0.
        while True:
            octree_depth += 1.
            scale = dimension / edge_size / 2 ** octree_depth
            if 0.010 < scale < 0.990:
                break
            if octree_depth > 9:
                raise BFException(ob, "Object too large for its edge size, split in parts.")
        return octree_depth, scale

class BFTInstantMesh(bpy.types.Operator):
    bl_idname="ops.bft_im"
    bl_label="Run InstantMesh"
    bl_description="Retopology of selected Objects with InstantMesh"
    bl_options={'REGISTER', 'UNDO'}
    bl_region_type="WINDOW"

    operation = bpy.props.StringProperty(name="operation", default='autorun')

    def execute(self, context):
        # Check and init
        if not bpy.context.selected_objects:
            self.report({'ERROR'}, "No selected Objects")
            return {'CANCELLED'}
        bft_im_filepath = bpy.path.abspath(context.user_preferences.addons[__name__].preferences.bft_im_filepath)
        if bft_im_filepath == "":
            self.report({'ERROR'}, "Set InstantMeshes filepath in File > User Preferences > Add-ons")
            return {'CANCELLED'}
        obj_filepath = bpy.app.tempdir + bpy.context.active_object.name + "_im.obj"
        ob_parent = bpy.context.active_object.parent
        wm = context.window_manager
        bpy.ops.view3d.snap_cursor_to_selected()
        # Export OBJ file
        try:
            bpy.ops.export_scene.obj(
                filepath=obj_filepath,
                use_selection=True,
                use_materials=False,
                use_mesh_modifiers=True,
                use_triangles=True,
                use_normals=False,
                use_edges=False,
                use_uvs=False,
            ) 
        except Exception as err:
            w.cursor_modal_restore()
            print(err)
            self.report({'ERROR'}, "Error while exporting OBJ file, see console")
            return {'CANCELLED'}
        # Execute IM
        w = context.window_manager.windows[0]
        w.cursor_modal_set("WAIT")
        creation_time = os.path.getmtime(obj_filepath)
        try:
            if self.operation == 'open':
                print(bft_im_filepath,              # IM binary
                    obj_filepath,                 # input file
#                    "--scale", str(wm.bft_edge),
                    "--faces", str(1000),
                    )
                subprocess.call([
                    bft_im_filepath,              # IM binary
                    obj_filepath,                 # input file
#                    "--scale", str(wm.bft_edge),  # desired edge length
                    "--faces", str(wm.bft_faces),
                    ])
            elif self.operation == 'autorun':
                subprocess.call([
                    bft_im_filepath,              # IM binary
                    "--output", obj_filepath,     # output file
                    "--threads", "2",             # calc threads
                    "--smooth",  "0",             # smoothing steps
#                    "--scale", str(wm.bft_edge),  # desired edge length
                    "--faces", str(wm.bft_faces),
#                   "--rosy",    "6",             # orientation sym type
#                   "--posy",    "6",             # position sym type
#                   "--deterministic",            # deterministic algorithms
                    obj_filepath,                 # input file
                    ])
        except Exception as err:
            w.cursor_modal_restore()
            print(err)
            self.report({'ERROR'}, "Error in InstantMeshes, see console")
            return {'CANCELLED'}
        w.cursor_modal_restore()
        # Check changed and import
        if(os.path.getmtime(obj_filepath) != creation_time):
            try:
                bpy.ops.import_scene.obj(filepath=obj_filepath)
            except Exception as err:
                print(err)
                self.report({'ERROR'}, "Error importing results back, see console")
                return {'CANCELLED'}
        else:
            self.report({'INFO'}, "No changes from original Object")
        # Remove tmp file
        try:
            os.remove(obj_filepath)
        except Exception as err:
            print(err)
        # Set Object origin to 3d cursor
        try:
            bpy.ops.object.origin_set(type='ORIGIN_CURSOR')
        except Exception as err:
            print(err)
            self.report({'ERROR'}, "Error setting Object origin, see console")
            return {'CANCELLED'}
        # Copy original properties to imported Object, set parent
        bpy.ops.object.bf_props_to_sel_obs()
        for ob in bpy.context.selected_objects:
            ob.parent = ob_parent
        return {'FINISHED'}


class BFTTopologyPanel(bpy.types.Panel):
    bl_label = "Topology"
    bl_idname = "OBJECT_PT_BFT"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'TOOLS'
    bl_category = "FDS Tools"
    bl_context = "objectmode"

    def draw(self, context):
        wm = context.window_manager
        ob = context.object
        layout = self.layout
        col = layout.column(align=True)
        col.label(text="Retopo:")
        col.operator("ops.bft_bl_shrinkwrap", text="Shrinkwrap").operation = "autorun"
        col.prop(wm, 'bft_edge')
        col.operator("ops.bft_bl_remesh", text="Remesh").operation = "autorun"
        col.prop(wm, 'bft_faces')
        col.operator("ops.bft_im", text="Run IM").operation = "autorun"
        col.operator("ops.bft_im", text="Open IM").operation = "open"
        col.operator("ops.bft_clean")
        

def register():
    bpy.utils.register_class(BFTClean)
    bpy.utils.register_class(BFTInstantMesh)
    bpy.utils.register_class(BFTBlShrinkwrap)
    bpy.utils.register_class(BFTBlRemesh)
    bpy.utils.register_class(BFTTopologyPanel)
    bpy.utils.register_class(BFTPref)
    
    bpy.types.WindowManager.bft_edge = FloatProperty(
        name='Desired Edge Length',
        description='Desired edge lenght',
        min=.000001, precision=2, step=100, default=1., unit='LENGTH',
        )
    bpy.types.WindowManager.bft_faces = IntProperty(
        name='Desired Faces',
        description='Desired number of faces',
        min=6, default=100,
        )

def unregister():
    bpy.utils.unregister_class(BFTClean)
    bpy.utils.unregister_class(BFTInstantMesh)
    bpy.utils.unregister_class(BFTBlShrinkwrap)
    bpy.utils.unregister_class(BFTBlRemesh)
    bpy.utils.unregister_class(BFTTopologyPanel)
    bpy.utils.unregister_class(BFTPref)
    
    try:
        del bpy.types.WindowManager.bft_edge
    except:
        pass


if __name__ == "__main__":
    register()
