def to_atom_instance(atoms):
    """
    Insted of creating a Model objct to store Atom objects, this function simply put Atom object into
    a dict. Use the following syntax to access the Atom object:
    Atom_dict[model][chain+seq.name], 
    e.g.: protein1[1][A12.CA], means access the Atom object in protein1, chain A, residue 12, atom CA.
    """
    atom_dict = {}               #Atom_dict: a dict contains dicts of atom object
    for i in range(len(atoms)):
        atom_dict[i+1] = {}
        for atom in atoms[i]:
            if atom[3] == ' ':
                atom[3] = 'A' #if chain id section is empty, then 'A' is arbitrarily assign as chain id.
                chain_id = atom[3]
                atom_name = atom[1]
                res_serial = str(atom[4])
            atom_dict[i+1][chain_id+res_serial+'.'+atom[1]]=Atom(atom) 
            #using chainid+ residue serial + atom name as identifier.
    return atom_dict