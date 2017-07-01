

def alignment_loader(alignment_file):
    """
    load the file and return the lines in this file.
    """
    seq_lines = []
    try:
        with open(alignment_file) as f:
            lines = f.readlines()
    except:
        print('Could not open pdb file!')
        raise
    for line in lines[1:]:
        if len(line) <= 1:
            pass
        else:
            seq_lines.append(line.strip('\n'))
    return seq_lines

def alignment_reader(alignment_file):
    '''
    format the sequence.
    {'header of seq 1':[sequences, 0, []],...}
    '''
    seq_lines = alignment_loader(alignment_file)
    line_len = len(seq_lines[0])
    seq_len = len(seq_lines[0].split()[1])
    header_len = line_len-seq_len  # to figure out how long if the header section

    aligned_seqs = {}
    for line in seq_lines:  #read each line, if they have the same header, then put them in the same item.
        if line is not None:
            header = line[:header_len]
            if header not in aligned_seqs:
                aligned_seqs[header] = []
                aligned_seqs[header].append(line[header_len:])
            else:
                aligned_seqs[header].append(line[header_len:])
    for s in aligned_seqs:
        aligned_seqs[s] = [''.join(aligned_seqs[s])]
    return aligned_seqs

def alignment_mapping(aligned_seqs):
    for s in aligned_seqs:
        aligned_seqs[s].append(0)  # for counting residues serial
        aligned_seqs[s].append([]) # for store the aligned residue serials
    seqs = sorted(i for i in aligned_seqs)
    alignment = aligned_seqs[seqs[0]][0]
    for i in range(len(alignment)):
        for s in seqs[1:]:
            seq = aligned_seqs[s][0][i]              # the single residue of the sequences
#            res_serial = aligned_seqs[s][1]        # the residue serial of the sequence
            align_seq_serial = aligned_seqs[s][2] # the list contains residue serial of the aligned part of sequence
            if seq !='-':
                aligned_seqs[s][1]+=1
                res_serial = aligned_seqs[s][1]
            if alignment[i] == '*':
                align_seq_serial.append(res_serial)
    del aligned_seqs[seqs[0]]
    return aligned_seqs

