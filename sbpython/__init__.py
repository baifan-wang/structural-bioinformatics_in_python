import numpy as np
from .Atom import *
from .Molecule import *
from .Residue import *
from .Model import *
from .pdb_paser import *
from .Geometry import *
from .Interaction import *
from .Superimposer import *
from .Structural_analysis import *
from .pyls import *
#from .atoms_to_model import *

def atoms_to_model(atoms_list):
    m = Model()
    i=0
    if len(atoms_list)==0:
        raise ValueError('No coordinates are found')
    for atoms in atoms_list:
        m1 = Molecule()
        i+=1
        m1.name = 'm'+str(i)
        r1 = Residue(atoms[0][4],atoms[0][2],atoms[0][3])
        for a in atoms:
            a1 = Atom(a)
            if a1.res_serial is r1.serial and a1.chain_id is r1.chain_id:
                r1.add_atom(a1)
            else:
                m1.add_residue(r1)
                r1 = Residue(a[4], a[2], a[3])
                r1.add_atom(a1)
        m1.add_residue(r1)
        m.add_molecule(m1)
    return m