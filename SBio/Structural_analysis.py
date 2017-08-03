from .Geometry import *
from .Molecule import *
import numpy as np

def get_residues(model, chain_id, residue_range):
    residues = []
    for c in chain_id:
        for i in range(residue_range[0],residue_range[1]+1):
            residue_id = c+str(i)
            residue = model.get_residue_by_id(residue_id)
            if residue is not None:
                residues.append(residue)
    return residues



def nucleic_acid_torsion(model, chain_id, residue_range):
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
    r = get_residues(model, chain_id, residue_range)
    try:
        for i in range(len(r)):
            if i-1 ==-1:
                alpha = None
                beta = None
            elif r[i-1].res_serial+1 != r[i].res_serial:
                alpha = None
                alpha = None
            else:
                alpha = get_torsion(r[i-1].O3_, r[i].P, r[i].O5_, r[i].C5_)
                beta = get_torsion(r[i].P, r[i].O5_, r[i].C5_, r[i].C4_)
            gamma = get_torsion(r[i].O5_, r[i].C5_, r[i].C4_, r[i].C3_)
            delta = get_torsion(r[i].C5_, r[i].C4_, r[i].C3_, r[i].O3_)
            if i+1 == len(r):
                epsilon = None
                zeta = None
            elif r[i].res_serial+1 != r[i+1].res_serial:  #deal with missing residues
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
            torsion.append([r[i].res_serial, alpha, beta, gamma, delta, epsilon, zeta, chi])
    except KeyError:
        pass
    return torsion

def nucleic_acid_torsion_plot(molecule, chain_id, residue_list):
    torsion = []
    for m in molecule.get_model():
        t = nucleic_acid_torsion(m,chain_id,residue_list)
        t = np.array(t)[:,1:]     #remove the first column, which is residue id.
        torsion.append(t)
    t1 = torsion[0]

    for t in torsion[1:]:
        t1 = np.vstack((t1,t))
    t1 = np.transpose(t1)
    return t1

def nucleic_acid_pucker(model, chain_id, residue_range):
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
    ps ={"C3'-endo":[0,36], #puckering_scheme
        "C4'-exo":[36,72],
        "C4'-endo":[72,108],
        "C1'-exo":[108,144],
        "C2'-endo":[144,180],
        "C3'-exo":[180,216],
        "C4'-endo":[216,252],
        "O4'-exo":[252,288],
        "C1'-endo":[288,324],
        "C2'-exo":[324,360]}
    r = get_residues(model, chain_id, residue_range)
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
                break
            pucker.append([i.res_serial, v0,v1,v2,v3,v4,tm,P,x])
    return pucker

def protein_tosion(model, chain_id, residue_range):
    """
    φ (phi):  C(i-1),N(i),CA(i),C(i) 
    ψ (psi):  N(i),CA(i),C(i),N(i+1) 
    """
    r = get_residues(model, chain_id, residue_range)
    torsion = []
    for i in range(len(r)):
        if i-1 == -1:
            phi = None
        elif r[i-1].res_serial+1 != r[i].res_serial:  #deal with missing residues
            phi = None
        else:
            phi = get_torsion(r[i-1].C, r[i].N, r[i].CA, r[i].C)
        if i+1 == len(r):
            psi = None
        elif r[i].res_serial+1 != r[i+1].res_serial:  #deal with missing residues
            phi = None
        else:
            psi = get_torsion(r[i].N, r[i].CA, r[i].C, r[i+1].N)
        torsion.append([r[i].res_serial, phi, psi])
    return torsion


def get_NOE(M1, M2, key='all'):
    peak_intense = {'all':[0,6.5], 'weak':[3.5,6.5],'medium':[2.6,5.0], 'strong':[1.8,3.6]}
    lcutoff = peak_intense[key][0]
    hcutoff = peak_intense[key][1]
    NOE = []
    for a1 in M1.get_atom():
        if a1.element == 'H':
            for a2 in M2.get_atom():
                if a2.element == 'H':
                    distance = get_distance(a1, a2)
                    if lcutoff<= distance <=hcutoff:
                        NOE.append([a1.id,a2.id,distance])
    return NOE