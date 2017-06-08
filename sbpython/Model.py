class Model():
    """
    Creat Model object.
    The Model object stores Molecule objects, i.e., the conformation(s) of a molecule or complex.
    Moleculeobjects are stored in a dict named molecules.
    """
    def __init__(self, molecules=None):
        if molecules == None:
            self._molecules = {}
        else:
            self._molecules = molecules
        self.id='model'

    def __str__(self):
        return self.id

    __repr__=__str__

    def __len__(self):
        """
        Return the number of conformations in this molecule.
        """
        return len(self._molecules)

    def __getattr__(self, args):
        """
        Enable the access Atom, Residue and Molecule in Model object with the following syntax:
        Model.Molecule.chian+Residue_sreial.Atom_name, e.g.: protein1.M1.A12.CA, protein1.M2.
        """
        if args in self._molecules: 
            return self._molecules[args]
        else:
            raise AttributeError('Not such molecule in this model!')

    def __contains__(self, molecule):
        """
        True if there is a molecule with the given id.
        """
        return (molecule.id in self._molecules) and (molecule.container is self)

    def add_molecule(self, molecule):
        """
        Adding a molecule to this molecule.
        raise KeyError if key conflicts.
        """
        id = molecule.id
        if id in self._molecules:
            raise KeyError("%s is already in this molecule!" %id)
        self._molecules[id] = molecule
        self._molecules[id].container = self

    def get_molecules(self):
        """
        return a generator of all molecules.
        """
        for m in self._molecules:
            yield self._molecules[m]



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