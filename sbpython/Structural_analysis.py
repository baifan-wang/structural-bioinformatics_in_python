from .Geometry import *
from .Molecule import *



def nucleic_acid_torsion(molecule, residue_range):
    """ 
        α, O3′−P−O5′−C5′
        β, P−O5′−C5′−C4′
        γ, O5′−C5′−C4′−C3′
        δ, C5′−C4′−C3′−O3′ 
        ε, C4′−C3′−O3′−P 
        ξ, C3′−O3′−P−O5′
        χ, O4′−C1′−N9−C4
        anti: −60°<=χ<=−180°,170°<=χ<=−180°
        syn:  30°<=χ<=90°
        """
    pass

def nucleic_acid_pucker(molecule, residue_range):
    """ 
        ν0, C4′−O4′−C1′−C2′
        ν1, O4′−C1′−C2′−C3′
        ν2, C1′−C2′−C3′−C4′
        ν4, C2′−C3′−C4′−O4′
        ν4, C3′−C4′−O4′−C1′
        tan(P) = (ν4+ν1-ν4-ν0)/2ν2(sin36°+sin72°)
        NORTH: 270° ≤ P ≤ 360°, 0° ≤ P < 90°
        SOUTH: 90° ≤ P < 270°
        """
    pass

def protein_tosion(molecule, residue_range):
    pass