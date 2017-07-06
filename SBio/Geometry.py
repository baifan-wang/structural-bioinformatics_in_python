import numpy as np

def vector_prod(a, b):
    """
    Compute vector product.
    """
    assert a.shape == (3, )
    assert b.shape == (3, )
    res = np.zeros(3, 'f')
    res[0] = a[1]*b[2] - a[2]*b[1]
    res[1] = a[2]*b[0] - a[0]*b[2]
    res[2] = a[0]*b[1] - a[1]*b[0]
    return res


def get_distance(c1, c2):
    """
    Compute the distance between point c1 and c2
    """
    c1 = c1.coord
    c2 = c2.coord
    assert c1.shape == (3, )
    assert c2.shape == (3, )
    d = c2 - c1
    return np.sqrt(np.sum(d*d))

def get_angle(c1, c2, c3):
    """
    Compute the angel between point c1, c2 and c3
    """
    c1 = c1.coord
    c2 = c2.coord
    c3 = c3.coord
    assert c1.shape == (3, )
    assert c2.shape == (3, )
    assert c3.shape == (3, )
    v1 = c1 - c2
    d1 = np.sqrt(np.sum(v1*v1))
    v2 = c3 - c2
    d2 = np.sqrt(np.sum(v2*v2))
    sca = np.dot(v1, v2)/(d1*d2)
    if sca<-1.0: 
        sca = -1.0
    elif sca>1.0: 
        sca = 1.0
    return np.arccos(sca)*180/np.pi

def vector_angle(a, b):
    d1 =np.sqrt(a.dot(a))
    d2 = np.sqrt(b.dot(b))
    sca = np.dot(a, b)/(d1*d2)

    if sca < -1.0: 
        sca = -1.0
    elif sca > 1.0: 
        sca = 1.0
    return np.arccos(sca)*180/np.pi

def get_torsion(c1, c2, c3, c4):
    """
    Compute the torsion angle between c1, c2, c3, c4.
    All coordinates are cartesian; result is in degrees.
    Raises a ValueError if angle is not defined.
    """

    c1 = c1.coord
    c2 = c2.coord
    c3 = c3.coord
    c4 = c4.coord
    assert c1.shape == (3, )
    assert c2.shape == (3, )
    assert c3.shape == (3, )
    assert c4.shape == (3, )
    tang=0.0
    a = c1-c2
    b = c3-c2
    c = vector_prod(a, b)

    a = c2-c3
    b = c4-c3
    d = vector_prod(a, b)

    dd = np.sqrt(np.sum(c*c))
    de = np.sqrt(np.sum(d*d))

    if dd<0.001 or de<0.001:
        raise ValueError ('Torsion angle undefined, degenerate points')

    vv = np.dot(c, d) / (dd*de);
    if vv<1.0: 
        tang = vv
    else: 
        tang = 1.0
    if tang < -1.0: 
        tang = -1.0
    tang = np.arccos(tang)
    tang = tang*57.296

    b = vector_prod(c, d)
    if np.dot(a, b) > 0.0: tang = -tang
    return tang


def get_plane(p1, p2, p3):
    """ 
    a = (p2.y-p1.y)*(p3.z-p1.z)-(p2.z-p1.z)*(p3.y-p1.y)
    b = (p2.z-p1.z)*(p3.x-p1.x)-(p2.x-p1.x)*(p3.z-p1.z)
    c = (p2.x-p1.x)*(p3.y-p1.y)-(p2.y-p1.y)*(p3.x-p1.x)
    d = (0-(a*p1.x+b*p1.y+c*p1.z))
    the norm of a plane is (a,b,c)
    """
    p1 = p1.coord
    p2 = p2.coord
    p3 = p3.coord
    assert p1.shape == (3, )
    assert p2.shape == (3, )
    assert p3.shape == (3, )
    a = (p2[1]-p1[1])*(p3[2]-p1[2])-(p2[2]-p1[2])*(p3[1]-p1[1])
    b = (p2[2]-p1[2])*(p3[0]-p1[0])-(p2[0]-p1[0])*(p3[2]-p1[2])
    c = (p2[0]-p1[0])*(p3[1]-p1[1])-(p2[1]-p1[1])*(p3[0]-p1[0])
    d = (0-(a*p1[0]+b*p1[1]+c*p1[2]))
    return (a, b, c, d)
