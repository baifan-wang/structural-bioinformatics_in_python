import numpy as np
class Residue():
    """Create Residue object.
    The Residue object stores residue serial, name, chain id and Atom objects.
    Atom objects are stored in a dict named atoms.
    """
    def __init__(self, serial, name, chain_id, atoms=None):
        self.container = 'M'            #indicate which molecule this atom belong to, object
        self.serial = serial            #residue serial number, int
        self.name = name                #residue name, str, eg., 'AlA'
        self.chain_id = chain_id        #chain id, str, single letter,
        if atoms == None:               #atoms should be a dict, deafult is None
            self.atoms = {}
        else:
            self.atoms = atoms 
        self.id = self.chain_id+str(self.serial)  #an id to identificate this residue

    def __str__(self):
        return ('Residue object: %s' %self.id)

    __repr__=__str__

    def __eq__(self, other):
        """
        Check whether this atom and other atom are the same
        """
        return (self.container == other.container) and (self.id == other.id) and (self.coordinates() == self.coordinates())

    def __len__(self):
        """
        return the number of atoms in this residue.
        """
        return len(self.atoms)

    def __contains__(self, atom):
        """
        True if there is a atom in this reidue.
        """
        return (atom in self.atoms.values())

    def get_atom(self):
        for a in self.atoms:
            yield self.atoms[a]

    def add(self, atom):
        """
        Add an Atom object.
        Checks for adding duplicate atoms, and raises a ValueError if so.
        """
        atname = atom.name
        atname = atname.replace("'", "_")  #use the '_' to replace " ' "
        if atname in self:
            raise ValueError("Atom %s defined twice in residue %s" % (atname, self))
        self.atoms[atname] = atom

    def __getattr__(self, args):
        if args in self.atoms:
            return self.atoms[args]
        else:
            raise AttributeError('Not such atom in this residue!')

    def atom_list(self):
        a = sorted(self.atoms[i].serial for i in self.atoms)
        return a

    def coordinates(self):
        """
        return the coordinates of this residue.
        """
        c= [i[0] for i in sorted([(self.atoms[i].coord, self.atoms[i].serial) for i in self.atoms],key=lambda x:x[1])]
        #looks ugly, there should be a better way
        #c = [(self.atoms[i].coord, self.atoms[i].serial) for i in self.atoms]
        #c = sorted(c, key=lambda x:x[1])
        #c = [i[0] for i in c]
        return np.array(c)
