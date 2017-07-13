import numpy as np
from .Structural_alignment import *
class Molecule:
    """
    Creat Molecule object.
    The Molecule object stores model objects, i.e., the conformation(s) of a model or complex.
    Moleculeobjects are stored in a dict named models.
    """
    def __init__(self, models=None):
        if models == None:
            self._models = {}
        else:
            self._models = models
            for m in self._models.vlaues():
                m._container = self
        self.name=''

    def __str__(self):
        if self.name is '':
            return "Molecule object: "
        else:
            return self.name

    __repr__=__str__

    def __len__(self):
        """
        Return the number of conformations in this model.
        """
        return len(self._models)

    def __getattr__(self, args):
        """
        Enable the access Atom, Residue and Molecule in Model object with the following syntax:
        Model.Molecule.chian+Residue_sreial.Atom_name, e.g.: protein1.M1.A12.CA, protein1.M2.
        """
        if args in self._models: 
            return self._models[args]
        else:
            raise AttributeError('Not such model in this model!')

    def __contains__(self, model):
        """
        True if there is a model with the given id.
        """
        return (model.name in self._models) and (model._container is self)

    def add_model(self, model):
        """
        Adding a model to this model.
        raise KeyError if key conflicts.
        """
        m_id = model.name
        if m_id in self._models:
            raise KeyError("%s is already in this model!" %m_id)
        self._models[m_id] = model
        self._models[m_id]._container = self

    def get_model(self):
        """
        return a generator of all models.
        """
        for m in self._models:
            yield self._models[m]

    def write_pdb(self, pdb):
        with open(pdb, 'w') as f:
            i = 1
            for model in self.get_model():
                f.write('MODEL       '+str(i)+'\n')
                for atom in model.get_atom():
                    f.write(atom.to_pdb()+'\n')
                f.write('TER     \n')
                f.write('ENDMDL\n')
                i+=1

    def getAverageCoords(self):
        """
        Compute the average coords of the alternate locations of an atom.
        """
        pass
