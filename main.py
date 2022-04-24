class FastaObj:
    def __init__(self, fasta_id, sequence):
        self.fasta_id = fasta_id
        self.sequence = sequence

    def get_id(self):
        return self.fasta_id

    def get_seq(self):
        return self.sequence

    def get_length(self):
        return len(self.sequence)


def create_fasta_list(file_name):
    with open(file_name) as file:
        current_id = ""
        current_seq = ""
        fasta_list = []

        for line in file:
            if line[0] == ">" and current_id == "":
                current_id = line[1:].rstrip()
            elif line[0] == ">" and current_id != "":
                fasta_list.append(FastaObj(current_id, current_seq))
                current_id = line[1:].rstrip()
                current_seq = ""
            else:
                current_seq += line.rstrip()
        else:
            fasta_list.append(FastaObj(current_id, current_seq))

    return fasta_list


def create_contig(fasta_list):
    contigs_list = []
    for fasta in fasta_list:
        current_seq = fasta.sequence
        # Want to find all other sequences with areas that overlap with the sequence we are checking.
        sequence_contigs = []
        for check_fasta in fasta_list:
            if fasta.get_id != check_fasta.get_id:
                check_seq = check_fasta.sequence
                for val, base in enumerate(current_seq):
                    sub_seq = current_seq[val:]
                    current_contig = ''
                    for i, (sub_seq_base, check_seq_base) in enumerate(zip(sub_seq, check_seq)):
                        if sub_seq_base == check_seq_base:
                            current_contig += sub_seq_base
                            if i == len(sub_seq) - 1 and len(current_contig) > 1:
                                sequence_contigs.append(f'{current_seq[:-len(current_contig)]}{check_seq}')
                        else:
                            break

        contigs_list += sequence_contigs
        print(current_seq, sequence_contigs)
    return contigs_list


seqs = create_fasta_list('sample.txt')
print(create_contig(seqs))
