import numpy as np
from numpy.linalg import svd, det


class SVDSuperimposer():
    """
    Class to run SVD alignment.
    """
    def __init__(self, ref_coords, target_coords):
        '''
        tatget coordinates (target_coords) will be put on top of reference coordinates(ref_coords).
        ref_coords and target_coords should be an Nx3 arrays
        '''
        self.ref_coords = ref_coords
        self.target_coords = target_coords

        assert ref_coords.shape[1] == 3
        assert target_coords.shape[1] == 3
        assert ref_coords.shape == target_coords.shape

        self.n = ref_coords.shape[0]
        self.trans_coords = None
        self.rotation = None
        self.trans = None
        self.rmsd = None
        self.init_rmsd = None

    def _compute_rmsd(self, coords1, coords2):
        """
        Return rmsd (root mean square deviation) between coords1 and coords2.
        """
        diff = coords1 - coords2
        n = coords1.shape[0]
        return np.sqrt(sum(sum(diff * diff)) / n)

    def run(self):
        """
        Superimpose the coordinate sets.
        """
        tar_coords = self.target_coords
        ref_coords = self.ref_coords
        # center on centroid
        av1 = sum(tar_coords) / self.n
        av2 = sum(ref_coords) / self.n
        tar_coords = tar_coords - av1
        ref_coords = ref_coords - av2
        # correlation matrix
        a = np.dot(np.transpose(tar_coords), ref_coords)
        u, d, vt = svd(a)
        self.rotation = np.transpose(np.dot(np.transpose(vt), np.transpose(u)))
        # check if we have found a reflection
        if det(self.rotation) < 0:
            vt[2] = -vt[2]
            self.rotation = np.transpose(np.dot(np.transpose(vt), np.transpose(u)))
        self.trans = av2 - np.dot(av1, self.rotation)

    def get_transformed(self):
        """
        Get the transformed coordinate set.
        y on x = np.dot(y, rot) + tran
        """
        if self.rotation is None:
            raise Exception("Nothing superimposed yet.")
        if self.trans_coords is None:
            self.trans_coords = np.dot(self.target_coords, self.rotation) + self.trans
        return self.trans_coords

    def get_rotran(self):
        """
        Right multiplying rotation matrix and translation.
        """
        if self.rotation is None:
            raise Exception("Nothing superimposed yet.")
        return self.rotation, self.trans

    def get_init_rmsd(self):
        """
        Root mean square deviation of untransformed coordinates.
        """
        if self.target_coords is None:
            raise Exception("No coordinates set yet.")
        if self.init_rmsd is None:
            self.init_rmsd = self._compute_rms(self.target_coords, self.ref_coords)
        return self.init_rmsd

    def get_rmsd(self):
        """
        Root mean square deviation of superimposed coordinates.
        """
        if self.rmsd is None:
            trans_coords = self.get_transformed()
            self.rmsd = self._compute_rmsd(trans_coords, self.ref_coords)
        return self.rmsd


