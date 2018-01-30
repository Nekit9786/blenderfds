import math
import sys
from functools import reduce

# increase the max number of recursive calls
sys.setrecursionlimit(10000) # my default is 1000, increasing too much may cause a seg fault

EPSILON = 1e-07
EPSILON2 = 1e-05

class Vector(object):
    "A 3D vector."
    def __init__(self, *args):
        self.x, self.y, self.z = 0., 0., 0.
        if len(args) == 3: # Vector(1,2,3)
            self.x = args[0]
            self.y = args[1]
            self.z = args[2]
        elif len(args) == 1:  # Vector([1,2,3])
            a = args[0]
            if isinstance(a, dict):
                self.x = a.get('x', 0.0)
                self.y = a.get('y', 0.0)
                self.z = a.get('z', 0.0)
            elif a is not None and len(a) == 3:
                self.x = a[0]
                self.y = a[1]
                self.z = a[2]

    def __repr__(self):
        return '({0}, {1}, {2})'.format(self.x, self.y, self.z)
            
    def clone(self):
        return Vector(self.x, self.y, self.z)
        
    def negated(self):
        return Vector(-self.x, -self.y, -self.z)

    def __neg__(self):
        return self.negated()
    
    def plus(self, a):
        return Vector(self.x+a.x, self.y+a.y, self.z+a.z)

    def __add__(self, a):
        return self.plus(a)
    
    def minus(self, a):
        return Vector(self.x-a.x, self.y-a.y, self.z-a.z)

    def __sub__(self, a):
        return self.minus(a)
    
    def times(self, a):
        return Vector(self.x*a, self.y*a, self.z*a)

    def __mul__(self, a):
        return self.times(a)

    def dividedBy(self, a):
        return Vector(self.x/a, self.y/a, self.z/a)

    def __truediv__(self, a):
        return self.dividedBy(float(a))

    def __div__(self, a):
        return self.dividedBy(float(a))
    
    def dot(self, a):
        return self.x*a.x + self.y*a.y + self.z*a.z
    
    def lerp(self, a, t):
        "Linear interpolation from self to a"
        return self.plus(a.minus(self).times(t));
    
    def length(self):
        return math.sqrt(self.dot(self))
    
    def unit(self):
        return self.dividedBy(self.length())
        
    def cross(self, a):
        return Vector(
            self.y * a.z - self.z * a.y,
            self.z * a.x - self.x * a.z,
            self.x * a.y - self.y * a.x)
          
    def __getitem__(self, key):
        return (self.x, self.y, self.z)[key]

    def __setitem__(self, key, value):
        l = [self.x, self.y, self.z]
        l[key] = value
        self.x, self.y, self.z = l
            
    def __len__(self):
        return 3
    
    def __iter__(self):
        return iter((self.x, self.y, self.z))
            
    def __repr__(self):
        return 'Vector(%.2f, %.2f, %0.2f)' % (self.x, self.y, self.z) 

    def __eq__(self, other):
        return \
            math.isclose(self.x, other.x, rel_tol=EPSILON) and \
            math.isclose(self.y, other.y, rel_tol=EPSILON) and \
            math.isclose(self.z, other.z, rel_tol=EPSILON)

    def __hash__(self):
        return hash((self.x, self.y, self.z))

    def isZero(self):
        return abs(self.x) < EPSILON and abs(self.y) < EPSILON and abs(self.z) < EPSILON
        
class Vertex(object):
    "A 3D vertex of a polygon."
    def __init__(self, pos, normal=None):
        self.pos = Vector(pos)
        self.normal = Vector(normal)
    
    def clone(self):
        return Vertex(self.pos.clone(), self.normal.clone())
    
    def flip(self):
        self.normal = self.normal.negated()

    def interpolate(self, other, t):
        """
        Return a new vertex between this vertex and `other` by linearly
        interpolating all properties using a parameter of `t`.
        """
        return Vertex(self.pos.lerp(other.pos,t),self.normal.lerp(other.normal,t))

    def isCollinear(self, v1, v2):
        """Return true if self and other vertices v1 and v2, all lie on the same line.
        >>> Vertex((-1.00, -1.00, -1.00)).isCollinear(Vertex((-1.00, 1.00, -1.00)),Vertex((1.00, 1.00, -1.00)))
        False
        >>> Vertex((-1.00, -1.00, -1.00)).isCollinear(Vertex((.0, .0, .0)),Vertex((1.00, 1.00, 1.00)))
        True
        """
        a, b, c = v1.pos, self.pos, v2.pos
        return (a - c).cross(a - b).isZero()

    def isInBBox(self, v1, v2):
        "Return true if self is in the bounding box (v1, v2)."
        p, q, r = v1.pos, self.pos, v2.pos
        return (p.x <= q.x <= r.x or r.x <= q.x <= p.x) and (p.y <= q.y <= r.y or r.y <= q.y <= p.y) and (p.z <= q.z <= r.z or r.z <= q.z <= p.z) 
        
    def isOnEdge(self, v1, v2):
        """Return true if self is on the edge (v1, v2).
        >>> Vertex((-0.5,-0.5,-1.0)).isOnEdge(Vertex((-0.5,1.0,-1.0)),Vertex((-0.5,-1.0, -1.0)))
        True
        """
        return self.isInBBox(v1, v2) and self.isCollinear(v1, v2)
    
    def __repr__(self):
        return repr(self.pos)

    def __eq__(self, other):
        return self.pos == other.pos
    
    def __hash__(self):
        return hash(self.pos)

class Plane(object):
    "A plane in 3D space."

    def __init__(self, normal, w):
        self.normal = normal
        self.w = w # perpendicular distance of the plane from (0, 0, 0)
    
    @classmethod
    def fromPoints(cls, a, b, c):
        n = b.minus(a).cross(c.minus(a)).unit()
        return Plane(n, n.dot(a))

    def clone(self):
        return Plane(self.normal.clone(), self.w)
        
    def flip(self):
        self.normal = self.normal.negated()
        self.w = -self.w

    def __repr__(self):
        return 'normal: {0} w: {1}'.format(self.normal, self.w)
    
    def splitPolygon(self, polygon, coplanarFront, coplanarBack, front, back):
        """
        Split `polygon` by this plane if needed, then put the polygon or polygon
        fragments in the appropriate lists. Coplanar polygons go into either
        `coplanarFront` or `coplanarBack` depending on their orientation with
        respect to this plane. Polygons in front or in back of this plane go into
        either `front` or `back`
        """
        COPLANAR = 0 # all the vertices are within EPSILON distance from plane
        FRONT = 1    # all the vertices are in front of the plane
        BACK = 2     # all the vertices are at the back of the plane
        SPANNING = 3 # some vertices are in front, some in the back

        # Classify each point as well as the entire polygon into one of the above classes.
        polygonType = 0
        vertexLocs = []
        
        numVertices = len(polygon.vertices)
        for i in range(numVertices):
            # Calc the distance between the vertex and this plane
            t = self.normal.dot(polygon.vertices[i].pos) - self.w
            # Classify the vertex
            loc = -1 # Is this necessary? FIXME
            if t < -EPSILON2: loc = BACK
            elif t > EPSILON2: loc = FRONT
            else: loc = COPLANAR
            polygonType |= loc
            vertexLocs.append(loc)
    
        # Put the polygon in the correct list, splitting it when necessary.
        if   polygonType == COPLANAR:
            if self.normal.dot(polygon.plane.normal) > 0: coplanarFront.append(polygon)
            else: coplanarBack.append(polygon)
        elif polygonType == FRONT: front.append(polygon)
        elif polygonType == BACK: back.append(polygon)
        elif polygonType == SPANNING: # the polygon is spanning, so cut it
            f, b = [], [] # front vertices, back vertices
            for i in range(numVertices):
                j = (i+1) % numVertices
                ti = vertexLocs[i] # can be BACK, COPLANAR, FRONT
                tj = vertexLocs[j]
                vi = polygon.vertices[i]
                vj = polygon.vertices[j]
                # the edge is on one side or spanning?
                if ti != BACK: f.append(vi)
                if ti != FRONT:
                    if ti != BACK: b.append(vi.clone())
                    else: b.append(vi)
                if (ti | tj) == SPANNING: # one is BACK, the other is FRONT, so split the edge
                    # interpolation weight at the intersection point
                    t = (self.w - self.normal.dot(vi.pos)) / self.normal.dot(vj.pos.minus(vi.pos))
                    # intersection point on the plane
                    v = vi.interpolate(vj, t)
                    f.append(v)
                    b.append(v.clone())
                    # add this cut to the list of polygon cuts, for t-junctions detection
                    polygon.cuts.append(v) # FIXME
            # Add cut polygons
            if len(f) >= 3: front.append(Polygon(f, polygon.shared, polygon.parent_csg, polygon.cuts))
            else: raise Exception('Wasting vertices (f)')
            if len(b) >= 3: back.append(Polygon(b, polygon.shared, polygon.parent_csg, polygon.cuts))
            else: raise Exception('Wasting vertices (b)')
                     
class Polygon(object):
    "Represents a convex polygon. The vertices must be coplanar and form a convex loop."
    
    def __init__(self, vertices, shared=None, parent_csg=None, cuts=None): # FIXME
        self.vertices = vertices
        self.shared = shared # FIXME delete
        self.parent_csg = parent_csg # FIXME
        self.cuts = cuts or [] # FIXME
        # Get plane from vertices, check if collinear triplet
        numVertices = len(self.vertices)
        #self.plane = Plane.fromPoints(vertices[0].pos, vertices[1].pos, vertices[2].pos) # Check for not collinear
        self.plane = None
        for i in range(numVertices-2):
            v0 = self.vertices[i]
            v1 = self.vertices[i+1]
            v2 = self.vertices[i+2]
            if not v0.isCollinear(v1, v2):
                self.plane = Plane.fromPoints(v0.pos, v1.pos, v2.pos)
                break
        if not self.plane:
            raise Exception('Zero area polygon, vertices:',v0,v1,v2)
        self.check()
        
    def clone(self):
        vertices = list(map(lambda v: v.clone(), self.vertices))
        return Polygon(vertices, self.shared, self.parent_csg, self.cuts) # FIXME
                
    def flip(self):
        self.vertices.reverse()
        map(lambda v: v.flip(), self.vertices)
        self.plane.flip()

    def __repr__(self): # FIXME add parent_csg
        return reduce(lambda x,y: x+y,
                      ['Polygon(['] + [repr(v) + ', ' \
                                       for v in self.vertices] + ['])'], '')

    def check(self): # FIXME
        # Degenerate
        if len(self.vertices) < 3: raise Exception('Polygon vertices < 3')
        # Not flat
        w = self.plane.w
        n = self.plane.normal
        for v in self.vertices:
            t = n.dot(v.pos) - w # distance from vertex to plane 
            if abs(t) > EPSILON: raise Exception('Polygon is not flat')
        # Not convex
        numVertices = len(self.vertices)
        pos0 = self.vertices[numVertices-2].pos
        pos1 = self.vertices[numVertices-2].pos
        for i in range(numVertices):
            pos2 = self.vertices[i].pos
            edge0 = pos1-pos0
            edge1 = pos2-pos1
            if edge0.cross(edge1).dot(n) < -EPSILON: raise Exception('Polygon is not convex')
            pos0 = pos1
            pos1 = pos2
        # Check area 0 FIXME because he has no plane, done in init
        pass
            
class BSPNode(object):
    """
    Holds a node in a BSP tree. A BSP tree is built from a collection of polygons
    by picking a polygon to split along. That polygon (and all other coplanar
    polygons) are added directly to that node and the other polygons are added to
    the front and/or back subtrees. This is not a leafy BSP tree since there is
    no distinction between internal and leaf nodes.
    """
    def __init__(self, polygons=None):
        self.plane = None # Plane instance
        self.front = None # BSPNode
        self.back = None  # BSPNode
        self.polygons = []
        if polygons: self.build(polygons)
            
    def clone(self):
        node = BSPNode()
        if self.plane: node.plane = self.plane.clone()
        if self.front: node.front = self.front.clone()
        if self.back:  node.back = self.back.clone()
        node.polygons = list(map(lambda p: p.clone(), self.polygons))
        return node
        
    def invert(self):
        """ 
        Convert solid space to empty space and empty space to solid space.
        """
        for poly in self.polygons:
            poly.flip()
        self.plane.flip()
        if self.front: 
            self.front.invert()
        if self.back: 
            self.back.invert()
        temp = self.front
        self.front = self.back
        self.back = temp
        
    def clipPolygons(self, polygons):
        """ 
        Recursively remove all polygons in `polygons` that are inside this BSP tree.
        """
        if not self.plane: return polygons[:]

        front = []
        back = []
        
        for poly in polygons: self.plane.splitPolygon(poly, front, back, front, back)

        if self.front: front = self.front.clipPolygons(front) # recurse on branches, conserve those polygons

        if self.back: back = self.back.clipPolygons(back) # recurse on branches
        else: back = [] # but remove polygons that are back of the leaves

        front.extend(back) # send all non removed polygons, for recursion
        return front
        
    def clipTo(self, bsp):
        """ 
        Remove all polygons in this BSP tree that are inside the other BSP tree.
        """
        self.polygons = bsp.clipPolygons(self.polygons)
        if self.front: self.front.clipTo(bsp)
        if self.back: self.back.clipTo(bsp)

    def fixTJunctions(self):
        self.fixPolygons(self.getAllCuts())

    def solderJunction(self):
        cuts = self.getNonManifoldVertices()
        self.fixPolygons(cuts)
        print("Non manifold vertices:\n", self.getNonManifoldVertices())

    def fixPolygons(self, cuts): # FIXME improve for speed
        # Check all polygons for all cuts
        for cut in cuts:
            for poly in self.allPolygons():
                numVertices = len(poly.vertices)
                for i in range(numVertices):
                    # Get the edge
                    v0 = poly.vertices[i]
                    v1 = poly.vertices[(i+1) % numVertices]
                    if v0.pos == cut.pos or v1.pos == cut.pos:
                        break
                    # If cut is on that edge, add to poly
                    if cut.isOnEdge(v0, v1):
                        poly.vertices.insert((i+1) % numVertices, cut.clone())
                        break
        
    def getAllCuts(self):
        # Get all cuts
        allCuts = set()
        for poly in self.allPolygons():
            allCuts.update(poly.cuts)
            poly.cuts = [] # FIXME not sure!
        return list(allCuts)
        
    def getNonManifoldEdges(self): # FIXME
        # Get all edges
        allEdges = []
        nonManifoldEdges = []
        for poly in self.allPolygons():
            numVertices = len(poly.vertices)
            allEdges.extend([(poly.vertices[i], poly.vertices[(i+1) % numVertices]) for i in range(numVertices)])
        # Check existence of antiEdge
        while allEdges:
            edge = allEdges.pop()
            antiEdge = (edge[1], edge[0])
            try:
                allEdges.remove(antiEdge)
            except:
                nonManifoldEdges.append(edge)
        return nonManifoldEdges

    def getNonManifoldVertices(self): # FIXME
        nonManifoldEdges = self.getNonManifoldEdges()
        nonManifoldVertices = set()
        nonManifoldVertices.update([edge[0] for edge in nonManifoldEdges])
        nonManifoldVertices.update([edge[1] for edge in nonManifoldEdges])
        return list(nonManifoldVertices)
        
    def allPolygons(self):
        """
        Return a list of all polygons in this BSP tree.
        """
        polygons = self.polygons[:]
        if self.front: 
            polygons.extend(self.front.allPolygons())
        if self.back: 
            polygons.extend(self.back.allPolygons())
        return polygons

    # def reTessellate(self):
        # polygons = self.polygons[:]
        # while polygons:
            # poly = polygons.pop()
            # numVertices = len(poly.vertices)
            # antiEdges = set((v[i].pos, v[(i+1) % numVertices]) for i in range(numVertices))

            # for p in polygons:
                # p_edges = set((v[(i+1) % numVertices], v[i].pos) for i in range(numVertices))
                # common_edge = antiEdges & p_edges
                # if common_edge:
                    # polygons.remove(p)
                    # poly.vertices = 
                    # polygons.append()
                
        # self.polygons = polygons
        
        # if self.front:
            # self.front.reTessellate()
        # if self.back:
            # self.back.reTessellate()            
        
    def build(self, polygons):
        """
        Build a BSP tree out of `polygons`. When called on an existing tree, the
        new polygons are filtered down to the bottom of the tree and become new
        nodes there. Each set of polygons is partitioned using the first polygon
        (no heuristic is used to pick a good split).
        """
        if len(polygons) == 0: return
        # first polygon plane is used as partition plane, if not existing
        if not self.plane: self.plane = polygons[0].plane.clone()
        # add first polygon to this node
        self.polygons.append(polygons[0])
        # split all other polygons using the first polygon's plane
        front = []
        back = []
        for poly in polygons[1:]:
            # coplanar front and back polygons go into self.polygons
            self.plane.splitPolygon(poly, self.polygons, self.polygons, front, back)
        # recursively build the BSP tree
        if len(front) > 0:
            if not self.front: self.front = BSPNode()
            self.front.build(front)
        if len(back) > 0:
            if not self.back: self.back = BSPNode()
            self.back.build(back)
            
            
if __name__ == "__main__":
    import doctest
    doctest.testmod()