import numpy as np
class Residue:
    """
    Create Residue object.
    The Residue object stores residue serial, name, chain id and Atom objects.
    Atom objects are stored in a dict named atoms.
    """
    def __init__(self, res_serial, name, chain_id, atoms=None):
        self._container  = ''            #indicate which molecule this residue belong to, object
        self.res_serial = res_serial            #residue serial number, int
        self.name = name                #residue name, str, eg., 'AlA'
        self.chain_id = chain_id        #chain id, str, single letter,
        if atoms == None:               #atoms should be a dict, deafult is None
            self._atoms = {}
        else:
            self._atoms = atoms             #if a dict contains Atom objects is suplied,
            for a in self._atoms.values():  #set the container of Atom objects to self.
                a_container  = self
        self.id = self.chain_id+'.'+self.name+'.'+str(self.res_serial)  #an id to identificate this residue
        self.serial = ''                # a absolute int to identify this residue
    def __str__(self):
        if self._container is '':
            return("Residue object: "+self.id)
        else:
            return str(self._container)+'.'+self.id

    __repr__=__str__

    def __eq__(self, other):
        """
        Check whether this atom and other atom are the same
        """
        return (self._container  == other._container ) and (self.id == other.id) and (self.coordinates() == self.coordinates())

    def __len__(self):
        """
        return the number of atoms in this residue.
        """
        return len(self._atoms)

    def __contains__(self, atom):
        """
        True if there is a atom in this reidue.
        """
        return (atom in self._atoms.values())

    def _sorted_atom_list(self):
        return sorted([i for i in self._atoms], \
            key=lambda x:self._atoms[x].serial)

    def get_atom(self):
        for a in self._sorted_atom_list():
            yield self._atoms[a]

    def add_atom(self, atom):
        """
        Add an Atom object.
        Checks for adding duplicate atoms, and raises a ValueError if so.
        """
        atname = atom.name
        atname = atname.replace("'", "_")  #use the '_' to replace " ' "
        if atname in self._atoms:
            raise ValueError("Atom %s defined twice in residue %s" % (atname, self))
        self._atoms[atname] = atom
        self._atoms[atname]._container = self

    def __getattr__(self, args):
        if args in self._atoms:
            return self._atoms[args]
        else:
            raise AttributeError('Not such atom in this residue!')

    def atom_list(self):
        a = sorted(self._atoms[i].atom_serial for i in self._atoms)
        return a

    def coordinates(self):
        """
        return the coordinates of this residue.
        """
        c= [i[0] for i in sorted([(self._atoms[i].coord, self._atoms[i].atom_serial) for i in self._atoms],key=lambda x:x[1])]
        #looks ugly, there should be a better way
        #c = [(self._atoms[i].coord, self._atoms[i].serial) for i in self._atoms]
        #c = sorted(c, key=lambda x:x[1])
        #c = [i[0] for i in c]
        return np.array(c)

    def bb_coordinates(self, mask=('N','CA','C',"P","O5'","C5'","C4'","C3'","O3'")):
        for a in self._atoms.values():
            if a.name in mask:
                yield a.coord


    def next_residue(self):
        """
        return the serial number of next residue.
        """
        pass