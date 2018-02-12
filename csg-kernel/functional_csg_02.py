import array
import sys
from bsp import *

# increase the max number of recursive calls
sys.setrecursionlimit(10000)


if __name__ == "__main__":
    import doctest
    doctest.testmod()
    
    name = "cube"
    
    geometry[0] = from_STL(filename='./test/{0}/{0}_a.stl'.format(name))
    geometry[1] = from_STL(filename='./test/{0}/{0}_b.stl'.format(name))

    a = BSP(igeom=0, ifaces=get_ifaces(0))
    b = BSP(igeom=1, ifaces=get_ifaces(1))

    clip_to(a, b)  # remove everything in a inside b
    clip_to(b, a)  # remove everything in b inside a

    geometry[0] = get_new_geom_from_bsp(a)
    geometry[1] = get_new_geom_from_bsp(b)

    # Send to STL
    to_STL(0, filename='./test/{0}/{0}_a_clipped.stl'.format(name))
    to_STL(1, filename='./test/{0}/{0}_b_clipped.stl'.format(name))

    print(a)
    print(b)
