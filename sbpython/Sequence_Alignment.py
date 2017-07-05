class Seq_align_parser():

    def __init__(self, alignment_file):

        self.align_f = alignment_file
        self.align_seqs = {}
        self.align_res_mask = {}

    def file_loader(self):
        """
        load the file and return the sequence lines in this file.
        """
        seq_lines = []   # list of the lines of sequence and marker in clustaw file.
        try:
            with open(self.align_f) as f:
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

    def alignment_reader(self):
        '''
        format the sequence.
        {'header of seq 1':[sequences, 0, []],...}
        '''
        seq_lines = self.file_loader()
        line_len = len(seq_lines[0])
        seq_len = len(seq_lines[0].split()[1])
        header_len = line_len-seq_len  # to figure out how long the header section is.

        for line in seq_lines:  #read each line, if they have the same header, put them in the same item.
            if line is not None:
                header = line[:header_len]
                if header not in self.align_seqs:
                    self.align_seqs[header] = []
                    self.align_seqs[header].append(line[header_len:])
                else:
                    self.align_seqs[header].append(line[header_len:])
        for s in self.align_seqs:
            self.align_seqs[s] = [''.join(self.align_seqs[s])]

    def alignment_mapping(self):
        """
        get the residue serial of the conserved residues for each sequence
        in the sequence alignment file.
        """
        for s in self.align_seqs:
            self.align_seqs[s].append(0)  # for counting residues serial
            self.align_seqs[s].append([]) # for store the aligned residue serials
        seqs = sorted(i for i in self.align_seqs)
        conserved_mask = self.align_seqs[seqs[0]][0]  # the conserved_mask using '*', ':', '.' to indicate
        for i in range(len(conserved_mask)):          # the residue conserved score.
            for s in seqs[1:]:
                seq = self.align_seqs[s][0][i]        # the single residue of the sequences
    #            res_serial = align_seqs[s][1]        # the residue serial of the sequence
                align_seq_serial = self.align_seqs[s][2] # the list contains residue serial of the aligned part of sequence
                if seq !='-':
                    self.align_seqs[s][1]+=1
                    res_serial = self.align_seqs[s][1]
                if conserved_mask[i] in ('*',):       # can use less conserved residue
                    align_seq_serial.append(res_serial)
        del self.align_seqs[seqs[0]]
        for i in self.align_seqs:
            self.align_res_mask[i] = self.align_seqs[i][2]

    def run(self):
        self.alignment_reader()
        self.alignment_mapping()
