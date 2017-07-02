class Structural_alignment()
    def __init__(self, molecule, residue_range, backbone_mark = True):
        self.molecule = molecule
        self.res_range = residue_range
        self.bb_mark = backbone_mark
    