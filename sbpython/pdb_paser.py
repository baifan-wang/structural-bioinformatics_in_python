import numpy as np
"""
PDB format:
01234567890123456789012345678901234567890123456789012345678901234567890123456789
ATOM    967 HO3'  DG B  15       2.235   9.325  11.866  1.00  0.00           H  
HETATM  863  K     K A 101      12.085  -0.820 -14.546  1.00 30.86           K  

abbrev.   type  column       detail
atype    str    [0:6]     Record Type
number   int    [6:11]    Atom serial number.
name     str    [11:16]   Atom name.
altLoc   str    [16]      Alternate location indicator.
res      str    [17:21]   Residue name.
chain    str    [21]      Chain identifier.
seq      int    [22:26]   Residue sequence number.
icode    str    [26]      Code for insertion of residues.
blank    str    [27:30]    
x        float  [30:38]    Orthogonal coordinates for X in Angstroms.
y        float  [38:46]    Orthogonal coordinates for Y in Angstroms.
z        float  [46:54]    Orthogonal coordinates for Z in Angstroms.
occu     str    [54:60]    Occupancy.
bfactor  str    [60:66]    Temperature factor.
blank2   str    [66:76]    
ele      str    [76:78]    Element symbol, right-justified.
charge   str    [78:80]    Charge on the atom. 
"""
def file_loader(file):
    """
    load the file and return the lines in this file.
    """
    try:
        with open(file) as f:
            lines = f.readlines()
    except:
        print('Could not open pdb file!')
        raise
    return lines

def pdb_to_atoms(pdb):
    """
    Read a pdb, and return a list contains information for each atom in pdb.
    atoms: store all of lines contain coordinates in pdb.
    """
    atoms, model_num = [[]],0
    lines = file_loader(pdb)
    identifiers = set(['ATOM  ', 'HETATM', 'ANISOU'])
    for line in lines:
        identifier = line[:6]
        if identifier in identifiers:
            atoms[model_num].append([int(line[6:11]), line[11:16].strip(" "),  line[17:21].strip(" "), line[21],\
             int(line[22:26]), float(line[30:38]), float(line[38:46]), float(line[46:54]), line[76:78].strip(' ')])
            #atoms: atom_serial, atom_name, residue_name, chain_id, residue_serial, x, y, z, element_type 
            #          0              1        2            3           4           5  6  7   8
        elif identifier=='ENDMDL':
            model_num += 1
            atoms.append([]) #put atoms in different in sublist of atoms
        else:
            pass
    while atoms[-1] == []:
        atoms.pop()
    for i in atoms: #if pdb doesn't contain element tag, guess from atom name
        if i[-1] =='':
            i[-1]=i[1][0]
    return atoms




