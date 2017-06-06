class Model():
    """
    Creat Model object.
    The Model object stores Molecule objects, i.e., the conformation(s) of a molecule or complex.
    Moleculeobjects are stored in a dict named molecules.
    """
    def __init__(self, molecules=None):
        if molecules == None:
            self.molecules = {}
        else:
            self.molecules = molecules
        self.id=''

    def __str__(self):
        return self.id

    __repr__=__str__

    def __len__(self):
        """
        Return the number of conformations in this molecule.
        """
        return len(self.molecules)

    def __getattr__(self, args):
        """
        Enable the access Atom, Residue and Molecule in Model object with the following syntax:
        Model.Molecule.chian+Residue_sreial.Atom_name, e.g.: protein1.M1.A12.CA, protein1.M2.
        """
        if args in self.molecules: 
            return self.molecules[args]
        else:
            raise AttributeError('Not such molecule in this model!')

    def __contains__(self, molecule):
        """
        True if there is a molecule with the given id.
        """
        id = molecule.id
        return (id in self.molecules) and (molecule.container == self)

    def add_residue(self, molecule):
        """Adding a molecule to this molecule.
        raise KeyError if key conflict.
        """
        id = molecule.id
        if id in self.residues:
            raise KeyError("%s is already in this molecule!" %id)
        self.molecules[id] = molecule
        self.molecules[id].container = self

    def all_molecules(self):
        """
        return a generator of all molecules.
        """
        pass

    def getAverageCoords(self):
        """
        Compute the average coords of the alternate locations of an atom.
        """
        if self.alternate != []:
            # get the alternate coords.
            coords = [self.coords,] + self.alternate.coords
            # compute the average as a list of coords.
            avCoords = (sum(coords)/len(coords)).tolist()
        else:
            # only one coords.
            avCoords = self.coords
        return avCoords