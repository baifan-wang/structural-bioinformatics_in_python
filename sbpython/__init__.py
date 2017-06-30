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
from .Plot import *
from .Alignment import *
#from .atoms_to_model import *

def atoms_to_model(atoms_list):
    m = Model()
    i=0
    if len(atoms_list)==0:
        raise ValueError('No coordinates are found')
    for atoms in atoms_list:
        m1 = Molecule()         # creat a Molecule object for each atom list in atoms_list
        i+=1                    # increasing model number by 1
        m1.name = 'm'+str(i)    # creat a name for each molecule
        r1 = Residue(atoms[0][4],atoms[0][2],atoms[0][3])  #create first residue from first atom.
        for a in atoms:         # for each line in atoms, which is a list
            a1 = Atom(a)        # create an Atom object form a list
            if a1.res_serial == r1.res_serial and a1.chain_id == r1.chain_id:
                r1.add_atom(a1)
            else:
                m1.add_residue(r1)
                r1 = Residue(a[4], a[2], a[3])
                r1.add_atom(a1)
        m1.add_residue(r1)
        m.add_molecule(m1)
    return m