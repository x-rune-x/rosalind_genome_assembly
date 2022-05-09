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


def create_contig(input_sequences, original_seqs):
    contigs_list = []
    for val, sequence in enumerate(input_sequences):
        if all(original_seq in sequence for original_seq in original_seqs):
            return sequence
        else:
            # Want to find all other sequences with areas that overlap with the sequence we are checking.
            sequence_contigs = []
            for check_sequence in input_sequences:
                if input_sequences.index(check_sequence) != val:
                    for index, base in enumerate(check_sequence):
                        sub_seq = sequence[index:]
                        current_contig = ''
                        for i, (sub_seq_base, check_seq_base) in enumerate(zip(sub_seq, check_sequence)):
                            if sub_seq_base == check_seq_base:
                                current_contig += sub_seq_base
                                if i == len(sub_seq) - 1 and len(current_contig) > 1:
                                    sequence_contigs.append(f'{sequence[:-len(current_contig)]}{check_sequence}')
                            else:
                                break
            if len(sequence_contigs):
                # Make assumption that the shortest contig has the most overlap.
                contigs_list.append(min(sequence_contigs))
            else:
                contigs_list.append(sequence)
    print(contigs_list)
    return create_contig(contigs_list, original_seqs)


sequence_list = [x.get_seq() for x in create_fasta_list('rosalind_long.txt')]
genome = create_contig(sequence_list, sequence_list)
print(genome)

with open('result.txt', 'w') as writer:
    writer.write(genome)

