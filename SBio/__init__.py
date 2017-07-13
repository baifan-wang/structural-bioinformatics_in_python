import numpy as np
import os
from .Atom import *
from .Molecule import *
from .Residue import *
from .Model import *
from .Parser import *
from .Geometry import *
from .Interaction import *
from .Superimposer import *
from .Structural_analysis import *
from .pyls import *
from .Plot import *
from .Sequence_Alignment import *
from .Structural_alignment import *
#from .atoms_to_model import *

def atoms_to_molecule(atoms_list):
    model = Molecule()
    i=0
    if len(atoms_list)==0:
        raise ValueError('No coordinates are found')
    for atoms in atoms_list:
        m1 = Model()         # creat a model object for each atom list in atoms_list
        i+=1                    # increasing model number by 1
        m1.name = 'm'+str(i)    # creat a name for each molecule
        m1.serial = i
        m = 1
        n = 1
        r1 = Residue(atoms[0][4],atoms[0][2],atoms[0][3])  #create first residue from first atom.
        r1.serial = m
        for a in atoms:         # for each line in atoms, which is a list
            a1 = Atom(a)        # create an Atom object form a list
            a1.serial = n
            if a1.res_serial == r1.res_serial and a1.chain_id == r1.chain_id:
                r1.add_atom(a1)
            else:
                m1.add_residue(r1)
                r1 = Residue(a[4], a[2], a[3])
                m+=1
                r1.serial = m
                r1.add_atom(a1)
            n+=1
        m1.add_residue(r1)
        model.add_molecule(m1)
    return model

def create_moelcule(file):
    if not os.path.exists(file):
        raise FileNotFoundError("can not find file: %s" %file)
    filename = os.path.split(file)[-1]
    ex = os.path.splitext(filename)[-1][1:]
    if ex not in parsers:
        print("The parser for this type of file: %s hasn't been implemented!" %ex)
        raise
    parser = parsers[ex]
    model = atoms_to_molecule(parser(file))
    model.name = filename
    return model