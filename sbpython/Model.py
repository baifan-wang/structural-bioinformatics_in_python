class Model():
    def __init__(self, molecules=None):
        if molecules == None:
            self.molecules = {}
        else:
            self.molecules = molecules
            
    def __len__(self):
        return len(self.molecules)

    def __getattr__(self, args):
        if args in self.molecules: 
            return self.molecules[args]
        else:
            raise AttributeError('Not such molecule in this model!')

    def __contains__(self, id):
        """
        True if there is a molecule with the given id.
        """
        return (id in self.molecules)

    def all_molecules(self):
        """
        return a generator of all molecules.
        """
        pass