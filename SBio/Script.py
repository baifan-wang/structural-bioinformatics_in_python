def find_g_quartet(mol):
    '''find the residues comprise the G-quartet or G-triple, G-base pair based on 
    the distance bewteen O6 and N1 atom of guanine residues'''
    g_res = []
    for res in mol.get_residue():
        if res.name in ['DG','DG3','DG5','BGM']: #BGM is a bromide guanine base
            g_res.append(res)
    g_pair = []
    for r1 in g_res:
        for r2 in g_res:
            if r1 is not r2:
                d = get_distance(r1.O6, r2.N1)
                if d<=3.2:
                    r1_id = r1.chain_id+'_'+r1.name+'_'+str(r1.res_serial)
                    r2_id = r2.chain_id+'_'+r2.name+'_'+str(r2.res_serial)
                    temp = [r1_id, r2_id]
                    for p in g_pair:
                        if r1_id in p:
                            p.append(r2_id)
                            temp.pop()
                        elif r2_id in p:
                            p.append(r1_id)
                            temp.pop()
                    if len(temp)==2:
                        g_pair.append([r1_id, r2_id])
    g_quartet = list((list(set(p)) for p in g_pair))
    return g_quartet # a list as following:
    #[['A_DG_14', 'A_DG_1', 'A_DG_20', 'A_DG_8'], \
    # ['A_DG_15', 'A_DG_19', 'A_DG_2', 'A_BGM_7'], \
    # ['A_DG_21', 'A_DG_9', 'A_DG_13']]