# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 13:46:27 2018

@author: eng4
"""

#def get_joined_geom(igeom0, igeom1):  # FIXME
#    verts = geometry[igeom0].verts[:]
#    faces = geometry[igeom0].faces[:]
#    nverts0 = int(len(verts)/3)
#    verts.extend(geometry[igeom1].verts[:])
#    for ivert in geometry[igeom1].faces[:]:
#        new_ivert = ivert + nverts0
#        faces.append(new_ivert)
#    return Geom(verts, faces)

#def remove_tjunctions(igeom, ifaces):  # FIXME edge loops ready!
#    """
#    Remove t-junctions from ifaces edges
#    >>> geometry[0] = Geom([-1,0,0, 0,0,0, 1,0,0, 1,1,0, -1,1,0, -1,-1,0, 0,-1,0, 1,-1,0],\
#                [1,3,4, 1,2,3, 0,5,6, 6,2,0, 6,7,2, 0,1,4,])  # open
#    >>> remove_tjunctions(igeom=0, ifaces=get_ifaces(igeom=0))
#    [[3, 4, 0, 5, 6, 7, 2], [1, 2, 0]]
#    """
#    border_halfedges = get_border_halfedges(igeom, ifaces)
#    edge_loops = get_edge_loops(igeom, border_halfedges)
#    edge_to_spl_ivert = {}
#    for edge_loop in edge_loops:
#        # Test all iverts from edge_loop for collinearity
#        # with other edges of the same edge_loop
#        iverts = edge_loop[:]
#        for ivert in iverts:
#            len_el = len(edge_loop)
#            v = Vector(get_vert(igeom, ivert))
#            for i in range(len_el):
#                el_ivert0 = edge_loop[(i)]
#                el_ivert1 = edge_loop[(i+1) % len_el]
#                if ivert == el_ivert0 or ivert == el_ivert1:
#                    continue
#                # Check if ivert is between
#                v0 = Vector(get_vert(igeom, el_ivert0))
#                v1 = Vector(get_vert(igeom, el_ivert1))
#                if v.is_within(v0, v1) and v.is_collinear(v0, v1):
#                    edge = (el_ivert0, el_ivert1)
#                    edge_to_spl_ivert[edge] = ivert
#    # Split
#    
#
#            for edge in spl_edges:
#                iface = border_halfedges[edge]
#            for ivert in spl_iverts:
#                pass
#            # do the splittings all together, using the same iverts
#            # append faces,  previous









#    geom_union(0, 1, name='test/cube/cube')

#    bsp_a = BSP(igeom=0, ifaces=get_ifaces(0))
##    bsp_b = BSP(igeom=1, ifaces=get_ifaces(1))
#    print("Initial")
#    print(bsp_a)
##    print(bsp_b)
#
#    # Clip
#    clip_to(bsp_a, bsp_b)  # remove everything in a inside b
#    print("After first clipping")
#    print(bsp_a)
#    print(bsp_b)
#
#    clip_to(bsp_b, bsp_a)  # remove everything in b inside a
#    print("After second clipping")
#    print(bsp_a)
#    print(bsp_b)
#
#    # Create new geometry
#    geometry[2] = get_new_geom_from_bsp(bsp_a)
#    geometry[3] = get_new_geom_from_bsp(bsp_b)
#
#    # Send to STL
#    to_STL(2, filename='{}_a_clipped.stl'.format(name))
#    to_STL(3, filename='{}_b_clipped.stl'.format(name))

#    geom_union(0, 1, name='cube_union')
#    geometry[0] = from_STL(filename='icosphere_a.stl')
#    geometry[1] = from_STL(filename='icosphere_b.stl')
#    geom_union(0, 1, name='icosphere_union')


#    # Create BSPs, this splits some triangles in the geometry
#    a = build_bsp(igeom=0, ifaces=get_ifaces(0))
#    b = build_bsp(igeom=1, ifaces=get_ifaces(1))
#    # Clip, this splits more triangles in the geometry
#    clip_to(a, b)
#    clip_to(b, a)
#    # Create new geometry
#    geometry.append(get_new_geom_from_bsp(a))
#    geometry.append(get_new_geom_from_bsp(b))
#    # Send to STL
#    to_STL(igeom=2, filename='a_clipped.stl')
#    to_STL(igeom=3, filename='b_clipped.stl')





