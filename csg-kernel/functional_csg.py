# import math
import array

geometry = list()

class Geom():
    def __init__(self, verts=None, faces=None, filename=None):
        self.verts = array.array('f')
        self.faces = array.array('i')
        if verts and faces:
            if (len(verts) % 3) != 0: raise Exception('verts length should be 3xn')
            self.verts.extend(verts)
            if (len(faces) % 3) != 0: raise Exception('faces length should be 3xn')
            self.faces.extend(faces)
        elif filename:
            self.from_STL(filename)
    
    def get_face(self, iface):
        """
        Get iface face connectivity
        >>> g = Geom([-1,-1,1, 1,-1,1, 1,1,1, -1,1,1], [0,1,2, 0,2,3])
        >>> g.get_face(1)
        array('i', [0, 2, 3])
        """
        return self.faces[3*iface:3*iface+3]

    def get_nfaces(self):
        """
        Get the range for faces
        >>> g = Geom([-1,-1,1, 1,-1,1, 1,1,1, -1,1,1], [0,1,2, 0,2,3])
        >>> g.get_nfaces()
        range(0, 2)
        """
        return range(int(len(self.faces)/3))
       
    def get_xyz_vert(self, ivert):
        """
        Get ivert vertex coordinates
        >>> g = Geom([-1,-1,1, 1,-1,1, 1,1,1, -1,1,1], [0,1,2, 0,2,3])
        >>> g.get_xyz_vert(2)
        array('f', [1.0, 1.0, 1.0])
        """
        return self.verts[3*ivert:3*ivert+3]

    def get_nverts(self):
        """
        Get the range for verts
        >>> g = Geom([-1,-1,1, 1,-1,1, 1,1,1, -1,1,1], [0,1,2, 0,2,3])
        >>> g.get_nverts()
        range(0, 4)
        """
        return range(int(len(self.verts)/3))
        
    def to_STL(self, filename):
        """
        Write self to STL file
        >>> g = Geom([-1,-1,1, 1,-1,1, 1,1,1, -1,1,1], [0,1,2, 0,2,3])
        >>> g.to_STL('doctest.stl')
        """
        with open(filename, 'w') as f:
            f.write('solid name\n')
            for iface in self.get_nfaces():
                f.write('facet normal 0 0 0\n')
                f.write('    outer loop\n')
                for ivert in self.get_face(iface):
                    f.write('        vertex {v[0]:.6f} {v[1]:.6f} {v[2]:.6f}\n'.format(v=self.get_xyz_vert(ivert)))
                f.write('    endloop\n')
                f.write('endfacet\n')
            f.write('endsolid name\n')
            
    def from_STL(self, filename):
        """
        Import verts and faces from STL file
        """
        import numpy
        from stl import mesh
        mesh = mesh.Mesh.from_file(filename)
        for iface, p in enumerate(mesh.points):
            self.verts.extend(p)
            self.faces.extend((3*iface, 3*iface+1, 3*iface+2))
        
        
if __name__ == "__main__":
    print("ciao")
    
    import doctest
    doctest.testmod()