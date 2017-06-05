import numpy as np
class Molecule():
    def __init__(self, residues=None, name=None, molecule_type=None):
        if residues is None:
            self.residues = {}
        else:
            self.residues = residues
        self.name = name
        self.molecule_type = molecule_type #protien, DNA, RNA, HYB, small, ion

    def __getattr__(self, args):
        if args in self.residues:
            return self.residues[args]
        else:
            raise AttributeError('Not such residue in this molecule!')

    def __len__(self):
        """Return the number of residues in this molecule."""
        return len(self.residues)

    def __contains__(self, id):
        """True if there is a residue in this molecule."""
        return (id in self.residues)

    def get_atoms(self):
        for r in self.get_residues():
            for a in r.atoms:
                yield r.atoms[a]

    def get_residues(self):
        for r in self.residues:
            yield self.residues[r]

    def residue_list(self):
        r = sorted(i for i in self.residues)
        return r

    def coordinates(self):
        """return the coordinates of all of the atoms"""
        c = []
        for r in self.residue_list():
            res = self.residues[r]
            for a in res.coordinates():
                c.append(a)
        return np.array(c)

    def subset_atoms(self, residue, mask):
        """
        return a subset of atoms define by te mask, such as backbone.
        """
        pass

    def center(self):
        """return the average position of all atoms."""
        pass

    def mass_center(self):
        """return the center of mass of all atoms."""
        pass

    def get_sequence(self):
        """return the sequence of this molecule."""
        pass

    def backbone(self):
        """
        return coordinate of backbone atoms of this molecule."""
        coordinate = []
        protein_bb = set('N', 'CA', 'C')
        nucleic_acid_bb = set("P","O5'","C5'","C4'","C3'","O3'")
        if self.molecule_type == 'protein':  #backbone atom of protein is 'N', 'CA', 'C':
            for a in self.get_atoms:
                if a.name in protein_bb:
                    coordinate.append(a.coord)
        elif self.molecule_type == 'DNA' or self.molecule_type == 'RNA':
            for a in self.get_atoms:
                if a.name in protein_bb:
                    coordinate.append(a.coord)
        else:
            print("This molecule type: %s doesn't have a defined backbone type." %molecule_type)
            return None
        return np.array(coordinate)