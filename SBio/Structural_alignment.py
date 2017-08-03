from .Model import *
from .Superimposer import *
import numpy as np
class Structural_alignment:
    def __init__(self, molecule_list, residue_range, chain_id='A', \
                       backbone_mark=True, update_coord=True, 
                       reference=0, ref_range=None):
        self.mol_list = molecule_list  # a list molecules to be aligned.
        self.res_range = residue_range # the serials of the residues to be superimposed.
        # should be a list composed of sub-list contains residues serial
        # must has the same length to molecule list
        # e.g., [[1,2,3], list(range(1,3))*3]
        assert len(self.mol_list) == len(self.res_range)
        self.chain_id = chain_id
        self.bb_mark = backbone_mark # if bb_mark == True, only backbone atom 
                                     # will be used for alignment
        self.update_coord = update_coord
        self.rmsd = []
        self.reference = reference
        self.ref_range = ref_range

    def get_coords(self, molecule, res_range):
        coords = []
        for r in res_range:
            residue_id = self.chain_id+str(r)
            residue = molecule.get_residue_by_id(residue_id)
            if self.bb_mark ==True:
                for c in residue.bb_coordinates():
                    coords.append(c)
            else:
                for c in residue.coordinates():
                    coords.append(c)
        return np.array(coords)

    def superimpose(self, ref_coords, tar_coords):
        sup = Superimposer(ref_coords, tar_coords)
        sup.run()
        rmsd = sup.get_rmsd()
        rotation, translation = sup.get_rotran()
        return rmsd, rotation, translation

    def set_reference(self):
        reference = self.reference
        if isinstance(reference, Model):
            self.ref_coords =  self.get_coords(self.reference, self.ref_range)
        elif isinstance(reference, int):
            self.ref_coords =  self.get_coords(self.mol_list[reference], self.res_range[reference])

    def run(self):
        self.set_reference()
        for m, r in zip(self.mol_list, self.res_range):
            tar_coords = self.get_coords(m, r)
            assert self.ref_coords.shape == tar_coords.shape
            rmsd, rot, trans = self.superimpose(self.ref_coords, tar_coords)
            self.rmsd.append(rmsd)
            if self.update_coord == True:
                m.update_atomic_coodinates(rot, trans)


