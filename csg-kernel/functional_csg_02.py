import array
import sys
from bsp import *

# increase the max number of recursive calls
sys.setrecursionlimit(10000)


#def get_joined_geom(igeom0, igeom1):  # FIXME
#    verts = geometry[igeom0].verts[:]
#    faces = geometry[igeom0].faces[:]
#    nverts0 = int(len(verts)/3)
#    verts.extend(geometry[igeom1].verts[:])
#    for ivert in geometry[igeom1].faces[:]:
#        new_ivert = ivert + nverts0
#        faces.append(new_ivert)
#    return Geom(verts, faces)


def remove_tjunctions(igeom, ifaces):  # FIXME edge loops ready!
    """
    """
    border_halfedges = get_border_halfedges(igeom, ifaces)
    edge_loops = get_edge_loops(igeom, border_halfedges)
    for edge_loop in edge_loops:
       for ivert in edge_loop:
           pass
           # if ivert collinear with any other edge of the loop,
           # add to the list of splittings of that edge
           # get the iface from border_halfedges
       for edge in spl_edges:
           iface = border_halfedges[edge]
           for ivert in spl_iverts:
           # do the splittings all together, using the same iverts
           # append faces, parent previous

def get_new_geom_from_bsp(bsp):  # FIXME test
#    """
#    Return a new Geom according to bsp contents
#    >>> geometry[0] = Geom([-1.0, -1.0, -1.0, -1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0, 1.0, -1.0, 1.0, 1.0, 1.0, 1.0, 1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0, -1.0], [0, 1, 2, 2, 3, 0, 3, 2, 4, 4, 5, 3, 5, 4, 6, 6, 7, 5, 1, 0, 7, 7, 6, 1, 7, 0, 3, 3, 5, 7, 4, 2, 1, 1, 6, 4])  # A cube
#    >>> bsp = BSP(igeom=0, ifaces=get_ifaces(0))
#    >>> geometry[0] = get_new_geom_from_bsp(bsp)
#    >>> get_nverts(igeom=0), get_nfaces(igeom=0)
#    (8, 12)
#    """
    # Init
    igeom = bsp.igeom
    selected_ifaces = set(get_all_ifaces_from_bsp(bsp))
    iface_to_parent = get_iface_to_parent(igeom)
    iface_to_children = get_iface_to_children(igeom)


    print("len(selected_ifaces) before replacement of parents:",len(selected_ifaces))

    # Choose the best faces by gruping siblings into parent
    print("len(selected_ifaces) before:",len(selected_ifaces))
    done = False
    while not done:
        done = True
        for iface in selected_ifaces:
            parent = iface_to_parent.get(iface, None)  # could be non existant
            if parent is None:
                continue
            siblings = set(iface_to_children[parent])
            if siblings and siblings <= selected_ifaces:
                selected_ifaces -= siblings  # remove siblings
                selected_ifaces.add(parent)  # and add parent instead
                done = False  # not finished
                break  # set is updated, restart of the loop needed
    print("len(selected_ifaces) after replacement of parents:",len(selected_ifaces))

    print("nverts before:",get_nverts(igeom))
    merge_duplicated_verts(igeom)
    print("nverts after:",get_nverts(igeom))

#    # Get iverts that deserve a check, those of faces with a parent,
#    # not replaced by the previous step
#    selected_ifaces_w_parent = [iface for iface in selected_ifaces if iface_to_parent.get(iface, None)]
#    check_double_iverts = []
#    for iface in selected_ifaces_w_parent:
#        check_double_iverts.extend(get_face(igeom, iface))
#    # Compare those iverts with the existing ones
#    nverts = get_nverts(igeom)
#    for check_double_ivert in check_double_iverts:
#        # Get the vert coordinates
#        check_double_vert = get_vert(igeom, check_double_ivert)
#        # Compare those with any vert
#        for ivert in range(nverts):
#            if (Vector(get_vert(igeom, ivert))-Vector(check_double_vert)).isZero():
#                # merge, check all faces for that ivert
#                for selected_iface in selected_ifaces:
#                    selected_face = get_face(igeom, selected_iface)
#                    if selected_face[0] == check_double_ivert:
#                        selected_face[0] = ivert
#                        print("Fix ivert:", ivert, check_double_ivert )
#                        break
#                    elif selected_face[1] == check_double_ivert:
#                        selected_face[1] = ivert
#                        print("Fix ivert:", ivert, check_double_ivert )
#                        break
#                    elif selected_face[2] == check_double_ivert:
#                        selected_face[2] = ivert
#                        print("Fix ivert:", ivert, check_double_ivert )
#                        break

    print("len(selected_ifaces) before merging of triangles:",len(selected_ifaces))

    # Merge siblings FIXME put in previous routine
    # The problem are TJunctions.
    done = False
    while not done:
        done = True
        for iface in selected_ifaces:
            parent = iface_to_parent.get(iface, None)  # could be non existant
            if parent is None:  # it is an original face
                continue
            descendants = get_iface_descendants(igeom, parent)
            descendants.remove(iface)
            for descendant in descendants:
                if descendant in selected_ifaces:
                    merged_iface = merge_ifaces(igeom, iface0=iface, iface1=descendant)
#                    print("iface:", iface, "descendant:", descendant, "merged:", merged_iface)
                    if merged_iface:
                        # remove merging couple from selected_ifaces
                        selected_ifaces.remove(iface)
                        selected_ifaces.remove(descendant)
                        # insert new iface in selected_ifaces
                        selected_ifaces.add(merged_iface)
                        done = False
                        break
            if not done:
                break
    print("len(selected_ifaces) after merging of triangles:",len(selected_ifaces))

    # Solve local TJunctions FIXME or solve them globally? let's see

    # Get only selected iverts, avoid duplication
    selected_iverts = []
    for iface in selected_ifaces:
        face = get_face(igeom, iface)
        selected_iverts.extend((face[0], face[1], face[2]))
    selected_iverts = set(selected_iverts)
    # Build the ivert to new_ivert dict, and the corresponding new_verts
    ivert_to_ivert = {}
    new_verts = []
    new_ivert = 0
    for ivert in selected_iverts:
        vert = get_vert(igeom, ivert)
        new_verts.extend(vert)
        ivert_to_ivert[ivert] = new_ivert
        new_ivert += 1
    # Build the new_faces
    new_faces = []
    for iface in selected_ifaces:
        face = get_face(igeom, iface)
        new_faces.extend((
                ivert_to_ivert[face[0]],
                ivert_to_ivert[face[1]],
                ivert_to_ivert[face[2]],
                ))
    print("nverts after removal of duplicates:", len(new_verts)/3)
    return Geom(new_verts, new_faces)


if __name__ == "__main__":
    import doctest
    doctest.testmod()

    geometry[0] = from_STL(filename='./test/sphere/sphere_a.stl')
    geometry[1] = from_STL(filename='./test/sphere/sphere_b.stl')
    geom_union(0, 1, name='test/sphere/sphere')


#    # Get geometries
#    print("Test cut: 1,0")
#
#    name="sphere"
#
#    geometry = [None, None, None, None]
#    geometry[0] = from_STL(filename='{}_a.stl'.format(name))
#    geometry[1] = from_STL(filename='{}_b.stl'.format(name))
#    print(geometry[0])
#    print(geometry[1])
#
#    bsp_a = BSP(igeom=0, ifaces=get_ifaces(0))
#    bsp_b = BSP(igeom=1, ifaces=get_ifaces(1))
#    print("Initial")
#    print(bsp_a)
#    print(bsp_b)
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
