from .Geometry import *
from .Atom import *
import numpy as np 

polar_atoms = set(['O','S','P','N',\
                      'F', 'Cl','Br','I',\
                      'K','Na','Li','Pb','Fe'])
#common polar element type in biomolecules as well as halogens and ions
def get_hydrogen_bond(donor,acceptor,donor_H=None,dis_cutoff=3.0,ang_cutoff=135):
    """
    Detect whether there is hydrogen bond interaction between the given atoms.
    It will calculate the distance between donor and acceptor atom(non-hydrogen),
    if hydrogen of the donor is provided, then the angle between donor-donor_H-acceptor
    will be calculated.
    The criteria for formation of hydrogen bond is that the distance is within the cutoff 
    and angle is larger than cutoff.
    The default distance cutoff and angle cutoff are adopt form Amber14 manual.
    """
    if donor.element not in polar_atoms or acceptor.element not in polar_atoms:
        print("Warning! Non-polar atom(s) detected!")
    dis = distance(donor, acceptor)
    if donor_H is not None:
        ang = angle(donor,donor_H,acceptor)
    else:
        ang = None
    if ang is not None:
        if dis <=dis_cutoff and ang >=ang_cutoff:
            return True, dis, ang
        else:
            return False, dis, ang
    else:
        if dis <=dis_cutoof:
            return True, dis, ang
        else:
            return False, dis, ang



def get_polar_interaction(a1, a2, dis_cutoff=3.2):

    dis = distance(a1,a2)
    if dis <=dis_cutoff:
        return True, dis
    else:
        return False, dis

def get_average_coords(atom_group):
    average_coords = sum(a.coord for a in atom_group)/len(atom_group)
    a = Dummy('average',average_coords)
    return a

def pi_pi_interaction(agroup1, agroup2):
    if len(agroup1)<3 or len(agroup2)<3:
        raise ValueError("At least three atoms are required for this type calculation")
    average_coords1 = get_average_coords(agroup1)
    average_coords2 = get_average_coords(agroup2)
    dis = distance(average_coords1, average_coords2)
    v1 = plane(agroup1[0],agroup1[1],agroup1[2])[0:3]
    v2 = plane(agroup2[0],agroup2[1],agroup2[2])[0:3]
    v1 = np.array(v1)
    v2 = np.array(v2)
    ang = vector_angle(v1, v2)
    return dis, ang


