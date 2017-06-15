import numpy as np
class Molecule():
    """
    Create Molecule object.
    The Molecule object stores Residue objects.
    Residue objects are stored in a dict named residues.
    """
    def __init__(self, residues=None, name=None, molecule_type=None):
        self._container = ''            #indicate which model this molecule belong to, object
        if residues is None:           #residues should be a dict, deafult is None
            self._residues = {}         #if None, cerat a empty dict
        else:
            self._residues = residues
            for r in self._residues.vlaues():
                r._container = self
        self.name = 'm1'                #molecule name, str, eg., 'M1'
        self.molecule_type = molecule_type #protien, DNA, RNA, HYB(hybird), small, ion, str
#        self.id = ''
        
    def __str__(self):
        if self._container is '':
            return("Molecule object: "+self.name)
        else:
            return str(self._container)+'.'+self.name

    __repr__=__str__


    def __getattr__(self, args):
        """
        Enable the access Atom and Residue in this Molecule object with the following syntax:
        Molecule.chian+Residue_sreial.Atom_name, e.g.: M1.A12.CA, M2.B1.
        """
        if args in self._residues:
            return self._residues[args]
        else:
            raise AttributeError('Not such residue in this molecule!')

    def __len__(self):
        """
        Return the number of residues in this molecule.
        """
        return len(self._residues)

    def __contains__(self, residue):
        """
        True if there is a residue in this molecule.
        """
        id = residue.chain_id+str(residue.serial)
        return (id in self._residues) and (residue._container == self)

    def add_residue(self, residue):
        """Adding a residue to this molecule.
        raise KeyError if key conflict.
        """
        id = residue.chain_id+str(residue.serial)
        if id in self._residues:
            raise KeyError("%s is already in this molecule!" %id)
        self._residues[id] = residue
        self._residues[id]._container = self

    def get_atom(self):
        for r in self.get_residue():
            for a in r.get_atom():
                yield a

    def get_residue(self):
        for r in self._residues:
            yield self._residues[r]

    def residue_list(self):
        r = sorted(i for i in self._residues)
        return r

    def coordinates(self):
        """
        return the coordinates of all of the atoms
        """
        c = []
        for r in self.residue_list():
            res = self._residues[r]
            for a in res.coordinates():
                c.append(a)
        return np.array(c)

    def subset_atoms(self, residue, mask):
        """
        return a subset of atoms define by te mask, such as backbone.
        """
        pass

    def center(self):
        """
        return the average position of all atoms.
        """
        pass

    def mass_center(self):
        """
        return the center of mass of all atoms.
        """
        pass

    def get_sequence(self):
        """
        return the sequence of this molecule.
        """
        pass

    def backbone(self):
        """
        return coordinate of backbone atoms of this molecule.
        """
        coordinate = []
        protein_bb = set(["N", "CA", "C", "O"])
        nucleic_acid_bb = set(["P","O5'","C5'","C4'","C3'","O3'"])
        if self.molecule_type == 'protein': 
            for a in self.get_atom():
                if a.name in protein_bb:
                    coordinate.append(a.coord)
        elif self.molecule_type == 'DNA' or self.molecule_type == 'RNA':
            for a in self.get_atom():
                if a.name in nucleic_acid_bb:
                    coordinate.append(a.coord)
        else:
            print("This molecule type: %s doesn't have a defined backbone atom." %molecule_type)
            return None
        return np.array(coordinate)

    def write_pdb(self, pdb):
        with open(pdb, 'w') as f:
            f.write('MODEL       1\n')
            for atom in self.get_atom():
                f.write(atom.to_pdb()+'\n')
            f.write('TER     \n')
            f.write('ENDMDL\n')