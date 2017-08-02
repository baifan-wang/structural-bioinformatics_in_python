import numpy as np
from .Atom import *
from .Data import *
class Model():
    """
    Create Model object.
    The Model object stores Residue objects.
    Residue objects are stored in a dict named residues.
    """
    def __init__(self, residues=None, name=None, molecule_type='protein'):
        self._container = ''            #indicate which molecule this model belong to, object
        if residues is None:           #residues should be a dict, deafult is None
            self._residues = {}         #if None, cerat a empty dict
        else:
            self._residues = residues
            for r in self._residues.vlaues():
                r._container = self
        self.name = 'm1'                #model name, str, eg., 'M1'
        self.molecule_type = molecule_type #protien, DNA, RNA, HYB(hybird), small, ion, str
#        self.id = ''
        self.serial = ''
    def __str__(self):
        if self._container is '':
            return("Model object: "+self.name)
        else:
            return str(self._container)+'.'+self.name

    __repr__=__str__


    def __getattr__(self, args):
        """
        Enable the access Atom and Residue in this Model object with the following syntax:
        Model.chian+Residue_sreial.Atom_name, e.g.: M1.A12.CA, M2.B1.
        """
        if args in self._residues:
            return self._residues[args]
        else:
            raise AttributeError('Not such residue in this molecule!')

    def __len__(self):
        """
        Return the number of residues and atoms in this molecule.
        """

        return len(self._residues)

    def get_atom_num(self):
        atom_num = 0
        for r in self._residues:
            atom_num+=len(self._residues[r])
        return atom_num

    def __contains__(self, residue):
        """
        True if there is a residue in this molecule.
        """
        r_id = residue.chain_id+str(residue.res_serial)
        return (r_id in self._residues) and (residue._container == self)

    def add_residue(self, residue):
        """
        Adding a residue to this model.
        raise KeyError if key conflict.
        """
        r_id = residue.chain_id+str(residue.res_serial)
        if r_id in self._residues:
            raise KeyError("%s is already in this model!" %r_id)
        self._residues[r_id] = residue
        self._residues[r_id]._container = self

    def _sorted_residue_list(self):
        return sorted([i for i in self._residues], \
            key=lambda x:self._residues[x].serial)

    def get_residue(self):
        for r in self._sorted_residue_list():
            yield self._residues[r]

    def get_atom(self):
        for r in self.get_residue():
            for a in r.get_atom():
                yield a

    def get_residue_by_id(self, residue_id):
        if residue_id in self._residues:
            return self._residues[residue_id]
        else:
            pass

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

    def subset_atoms(self, chain='A',residue_resial=None, residue_name=None, atom_serial=None, atom_name=None):
        """
        return a subset of atoms define by te mask, such as backbone.
        """
        pass

    def center(self):
        """
        return the average position of all atoms and pack into a Dummt atom object.
        """
        center_coord = sum(self.coordinates())/len(self.coordinates())
        a = Dummy('Center of '+self.name, center_coord)
        return a

    def mass_center(self):
        """
        return the center of mass of all atoms.
        """
        pass

    def get_sequence(self, chain_id):
        """
        return the sequence of this molecule.
        """
        sequence = []
        temp_residue = []
        for r in self._residues:
            if r[0] == chain_id:
                temp_residue.append(r)
        temp_residue = sorted(temp_residue, key = lambda x:int(x[1:]) )
        for i in temp_residue:
            residue = self._residues[i]
            if residue.name in aa_res:
                sequence.append(aa_res[residue.name])
        return ''.join(sequence)

    def write_fasta(self, chain_id, fasta_file, comment=None):
        if comment is None:
            comment = 'Sequence extract from %s' %self
        sequence = self.get_sequence(chain_id)
        header = '>'+self.name+':'+chain_id+"|PDBID|CHAIN|SEQUENCE " + comment
        with open(fasta_file, 'w') as f:
            f.write(header+'\n')
            i = 0
            l = 80
            while i<len(sequence):
                f.write(sequence[i:l]+'\n')
                i+=80
                l+=80
        print('Writing sequence to %s' %fasta_file)

    def get_bb_coords(self):
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
            print("This molecule type: %s doesn't have defined backbone atoms." %molecule_type)
            return None
        return np.array(coordinate)

    def write_pdb(self, pdb):
        with open(pdb, 'w') as f:
            f.write('MODEL       1\n')
            for atom in self.get_atom():
                res = atom.res_serial
                break   #get first residue serial number
            for atom in self.get_atom():
                if (atom.res_serial != res+1) and (atom.res_serial != res):
                    f.write('TER     \n')   
                f.write(atom.to_pdb()+'\n')
            f.write('TER     \n')
            f.write('ENDMDL\n')

    def update_atomic_coodinates(self, rotation, translation):
        """
        update the coordinates of each atom in this molecule
        with the rotation matrix and translation vector.
        """
        for a in self.get_atom():
            a.transform(rotation, translation)

    def mw(self):
      """
      Returns molecular weight
      """
      mw = 0
      for a in self.get_atom():
        mw += atomic_weight[a.element]
      return mw

    def get_dist_matrix(self, savefile=None):
        """
        Compute the complete inter-atomic distance matrix, which has order of n2 time complexity!
        """
        atom_num = self.get_atom_num()
        self.dist_matrix = np.zeros((atom_num, atom_num)) #creat a zeroed atom_num*atom_num array
        for a1 in self.get_atom():
            i = a1.serial-1
            for a2 in self.get_atom():
                j = a2.serial-1
                self.dist_matrix[i,j] = a1.distance(a2)
        if savefile is not None:
            with open('savefile', 'w') as f:
                for i in self.dist_matrix:
                    for a in i:
                        data = '%.3f' %a
                        f.write(str(data)+'\t')
                    f.write('\n')