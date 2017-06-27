from .Geometry import *
from .Molecule import *
import numpy as np

def get_residues(molecule, chain_id, residue_range):
    residues = []
    for c in chain_id:
        for i in range(residue_range[0],residue_range[1]+1):
            residue_id = c+str(i)
            residues.append(molecule.get_residue_by_id(residue_id))
    return residues


def nucleic_acid_torsion(molecule, chain_id, residue_range):
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
    torsion = []
    r = get_residues(molecule, chain_id, residue_range)
    try:
        for i in range(len(r)):
            if i-1 ==-1:
                alpha = None
                beta = None
            else:
                alpha = get_torsion(r[i-1].O3_, r[i].P, r[i].O5_, r[i].C5_)
                beta = get_torsion(r[i].P, r[i].O5_, r[i].C5_, r[i].C4_)
            gamma = get_torsion(r[i].O5_, r[i].C5_, r[i].C4_, r[i].C3_)
            delta = get_torsion(r[i].C5_, r[i].C4_, r[i].C3_, r[i].O3_)
            if i+1 == len(r):
                epsilon = None
                zeta = None
            else:
                epsilon = get_torsion(r[i].C4_, r[i].C3_, r[i].O3_, r[i+1].P)
                zeta = get_torsion(r[i].C3_, r[i].O3_, r[i+1].P, r[i+1].O5_)
            if r[i].name in set(['A','DA','ADE','G','DG','GUA']):#purine base
                chi = get_torsion(r[i].O4_,r[i].C1_,r[i].N9,r[i].C4)
            elif r[i].name in set(['T','DT','THY','C','DC','CYT']):#purine base
                chi = get_torsion(r[i].O4_,r[i].C1_,r[i].N1,r[i].C2)
            else:
                chi = None
            torsion.append([alpha, beta, gamma, delta, epsilon, zeta,chi])
    except KeyError:
        pass
    return torsion

def nucleic_acid_pucker(molecule, chain_id, residue_range):
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
    ps ={"C3'-endo":[0.0001,36], #puckering_scheme
        "C4'-exo":[36.0001,72],
        "C4'-endo":[72.0001,108],
        "C1'-exo":[108.0001,144],
        "C2'-endo":[144.0001,180],
        "C3'-exo":[180.0001,216],
        "C4'-endo":[216.0001,252],
        "O4'-exo":[252.0001,288],
        "C1'-endo":[288.0001,324],
        "C2'-exo":[324.0001,360]}
    r = get_residues(molecule, chain_id, residue_range)
    pucker = []
    Pconst = np.sin(np.pi/5) + np.sin(np.pi/2.5)
    for i in r:
        v0 = get_torsion(i.C4_, i.O4_, i.C1_, i.C2_)
        v1 = get_torsion(i.O4_, i.C1_, i.C2_, i.C3_)
        v2 = get_torsion(i.C1_, i.C2_, i.C3_, i.C4_)
        v3 = get_torsion(i.C2_, i.C3_, i.C4_, i.O4_)
        v4 = get_torsion(i.C3_, i.C4_, i.O4_, i.C1_)
        tan_p = (v4+v1-v3-v0)/2*v2*(np.sin(0.2)+np.sin(0.4))
        P0 = np.arctan2((v4 + v1 - v3 - v0),(2.0 * v2 * Pconst))
        tm = v2/np.cos(P0)  # amplitude
        P = np.degrees(P0)
        if P<0:
            P = P+360
        for x in ps:
            if ps[x][0] <=P <= ps[x][1]:
                pucker.append([v0,v1,v2,v3,v4,tm,P,x])
    return pucker

def protein_tosion(molecule, chain_id, residue_range):
    """
    φ (phi):  C(i-1),N(i),CA(i),C(i) 
    ψ (psi):  N(i),CA(i),C(i),N(i+1) 
    """
    r = get_residues(molecule, chain_id, residue_range)
    torsion = []
    for i in range(len(r)):
        if i-1 == -1:
            phi = None
        else:
            phi = get_torsion(r[i-1].C, r[i].N, r[i].CA, r[i].C)
        if i+1 == len(r):
            psi = None
        else:
            psi = get_torsion(r[i].N, r[i].CA, r[i].C, r[i+1].N)
        torsion.append([phi,psi])
    return torsion

def plot_phi_psi(phi_psi):
    torsion = np.array(phi_psi[1:-1])
    pass