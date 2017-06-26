rom .Geometry import *
from .Molecule import *
import numpy as np


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
    r = molecule.get_subset(residue_range)
    torsion = []
    try:
        for i in range(len(r)):
            if i-1 ==-1:
                alpha = None
            else:
                alpha = dihedral(r[i-1].O3_, r[i].P, r[i].O5_, r[i].C5_)
            beta = dihedral(r[i].P, r[i].O5_, r[i].C5_, r[i].C4_)
            gamma = dihedral(r[i].O5_, r[i].C5_, r[i].C4_, r[i].C3_)
            delta = dihedral(r[i].C5_, r[i].C4_, r[i].C3_, r[i].O3_)
            if i+1 == len(r):
                epsilon = None
                zeta = None
            else:
                epsilon = dihedral(r[i].C4_, r[i].C3_, r[i].O3_, r[i+1].P)
                zeta = dihedral(r[i].O3_, r[i].C4_, r[i].C3_, r[i+1].P, r[i+1].O3_)
            if r[i].name in set(['A','DA','ADE','G','DG','GUA']):#purine base
                chi = dihedral(r[i].O4_,r[i].C1_,r[i].N9,r[i].C4)
            elif r[i].name in set(['T','DT','THY','C','DC','CYT']):#purine base
                chi = dihedral(r[i].O4_,r[i].C1_,r[i].N1,r[i].C2)
            else:
                chi = None
            torsion.append([alpha, beta, gamma, delta, epsilon, zeta,chi])
    except KeyError:
        pass
    return torsion

def nucleic_acid_pucker(molecule, residue_range):
    """ 
        ν0, C4′−O4′−C1′−C2′
        ν1, O4′−C1′−C2′−C3′
        ν2, C1′−C2′−C3′−C4′
        ν43, C2′−C3′−C4′−O4′
        ν4, C3′−C4′−O4′−C1′
        tan(P) = (ν4+ν1-ν3-ν0)/2ν2(sin36°+sin72°), 36° = 0.2arc
        NORTH: 270° ≤ P ≤ 360°, 0° ≤ P < 90°
        SOUTH: 90° ≤ P < 270°
        """
    pucker = []
    Pconst = np.sin(np.pi/5) + np.sin(np.pi/2.5)
    for i in r:
        v0 = dihedral(i.C4_, i.O4_, i.C1_, i.C2_)
        v1 = dihedral(i.O4_, i.C1_, i.C2_, i.C3_)
        v2 = dihedral(i.C1_, i.C2_, i.C3_, i.C4_)
        v3 = dihedral(i.C2_, i.C3_, i.C4_, i,O4_)
        v4 = dihedral(i.C3_, i.C4_, i,O4_, i.C1_)
        tan_p = (v4+v1-v3-v0)/2*v2*(np.sin(0.2)+np.sin(0.4))
        P0 = np.arctan2((v4 + v1 - v3 - v0),(2.0 * v2 * Pconst))
        tm = v2/np.cos(P0)  # amplitude
        P = np.degrees(P0)
        if P<0:
            P = P+360
        p.append([v0,v1,v2,v3,v4,P])
    return pucker

def protein_tosion(molecule, residue_range):
    """
    φ (phi):  C(i-1),N(i),CA(i),C(i) 
    ψ (psi):  N(i),CA(i),C(i),N(i+1) 
    """
    torsion = []
    for i in range(len(r)):
        if i-1 == -1:
            phi = None
        else:
            phi = dihedral(r[i-1].C, r[i].N, r[i].CA, r[i].C)
        if i+1 == len(r):
            psi = None
        else:
            psi = dihedral(r[i].N, r[i].CA, r[i].C, r[i+1].N)
        torsion.append([phi,psi])

    return torsion