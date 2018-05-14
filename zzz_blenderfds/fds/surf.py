"""BlenderFDS, FDS SURF routines."""

import bpy

predefined = ("INERT", "OPEN", "MIRROR", "HVAC", "PERIODIC")


def has_predefined():
    """Check predefined materials."""
    return set(predefined) <= set(bpy.data.materials.keys())


def set_predefined(context):
    """Set BlenderFDS predefined materials/bcs."""
    mas = bpy.data.materials.keys()
    value = str()
    if "INERT" not in mas:
        value += "&SURF ID='INERT' RGB=204,204,51 /\n"
    if "OPEN" not in mas:
        value += "&SURF ID='OPEN' RGB=51,204,204 TRANSPARENCY=.2 /\n"
    if "PERIODIC" not in mas:
        value += "&SURF ID='PERIODIC' RGB=204,204,204 TRANSPARENCY=.2 /\n"
    if "MIRROR" not in mas:
        value += "&SURF ID='MIRROR' RGB=51,51,204 /\n"
    if "HVAC" not in mas:
        value += "&SURF ID='HVAC' RGB=51,51,204 /\n"
    if value:
        context.scene.from_fds(context, value)
