import numpy as np
class Residue():
    def __init__(self, serial, name, chain_id, atoms=None):
        self.serial = serial
        self.name = name
        self.chain_id = chain_id
        if atoms == None:
            self.atoms = {}
        else:
            self.atoms = atoms 
        self.id = self.chain_id+str(self.serial)

    def __str__(self):
        return ('Residue object: %s' %self.id)

    __repr__=__str__

    def __len__(self):
        """
        return the number of atoms in this residue.
        """
        return len(self.atoms)

    def __contains__(self, atom_name):
        """
        True if there is a atom in this reidue.
        """
        return (atom_name in self.atoms)

    def get_atom(self):
        for a in self.atoms:
            yield self.atoms[a]

    def add(self, atom):
        """
        Add an Atom object.
        Checks for adding duplicate atoms, and raises a ValueError if so.
        """
        atname = atom.name
        if atname.endswith("'"):      #some atoms have " ' " in they names.
            atname = atname[:-1]+'_'  #use the '_' to replace " ' "
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
