import numpy as np
from numpy.linalg import svd, det


class SVDSuperimposer():
    """Class to run SVD alignment,

    SVDSuperimposer finds the best rotation and translation to put
    two point sets on top of each other (minimizing the RMSD). 

    SVD stands for Singular Value Decomposition, which is used to calculate
    the superposition.

    Reference:

    Matrix computations, 2nd ed. Golub, G. & Van Loan, CF., The Johns
    Hopkins University Press, Baltimore, 1989
    """
    def __init__(self, ref_coords, coords):
        '''
            coords will be put on top of ref_coords.

            - ref_coords: an NxDIM array
            - coords: an NxDIM array

            DIM is the dimension of the points, N is the number
            of points to be superimposed.'''
        self.ref_coords = ref_coords
        self.coords = coords


        assert ref_coords.shape[1] == 3
        assert coords.shape[1] == 3
        assert ref_coords.shape == coords.shape

        self.n = ref_coords.shape[0]
        self.transformed_coords = None
        self.rot = None
        self.tran = None
        self.rms = None
        self.init_rms = None

    def _compute_rms(self, coords1, coords2):
        """
        Return rms deviations between coords1 and coords2.
        """
        diff = coords1 - coords2
        l = coords1.shape[0]
        return np.sqrt(sum(sum(diff * diff)) / l)

    def run(self):
        """Superimpose the coordinate sets."""
        coords = self.coords
        ref_coords = self.ref_coords
        # center on centroid
        av1 = sum(coords) / self.n
        av2 = sum(ref_coords) / self.n
        coords = coords - av1
        ref_coords = ref_coords - av2
        # correlation matrix
        a = np.dot(np.transpose(coords), ref_coords)
        u, d, vt = svd(a)
        self.rot = np.transpose(np.dot(np.transpose(vt), np.transpose(u)))
        # check if we have found a reflection
        if det(self.rot) < 0:
            vt[2] = -vt[2]
            self.rot = np.transpose(np.dot(np.transpose(vt), np.transpose(u)))
        self.tran = av2 - np.dot(av1, self.rot)

    def get_transformed(self):
        """
        Get the transformed coordinate set.
        """
        if self.coords is None or self.ref_coords is None:
            raise Exception("No coordinates set.")
        if self.rot is None:
            raise Exception("Nothing superimposed yet.")
        if self.transformed_coords is None:
            self.transformed_coords = np.dot(self.coords, self.rot) + self.tran
        return self.transformed_coords

    def get_rotran(self):
        """
        Right multiplying rotation matrix and translation.
        """
        if self.rot is None:
            raise Exception("Nothing superimposed yet.")
        return self.rot, self.tran

    def get_init_rms(self):
        """
        Root mean square deviation of untransformed coordinates.
        """
        if self.coords is None:
            raise Exception("No coordinates set yet.")
        if self.init_rms is None:
            self.init_rms = self._compute_rms(self.coords, self.ref_coords)
        return self.init_rms

    def get_rms(self):
        """
        Root mean square deviation of superimposed coordinates.
        """
        if self.rms is None:
            transformed_coords = self.get_transformed()
            self.rms = self._compute_rms(transformed_coords, self.ref_coords)
        return self.rms
