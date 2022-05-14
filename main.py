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
    new_contigs = []
    for val, sequence in enumerate(input_sequences):
        if all(original_seq in sequence for original_seq in original_seqs):
            return sequence
        else:
            # Want to find all other sequences with areas that overlap with the sequence we are checking.
            for check_sequence in input_sequences:
                # Prevent checking string against itself but don't exclude duplicate sequences.
                if input_sequences.index(check_sequence) != val:
                    combined_strings = combine_strings_with_overlap(sequence, check_sequence)
                    if combined_strings:
                        new_contigs.append(combined_strings)
                        input_sequences.remove(sequence)
                        input_sequences.remove(check_sequence)

    if len(new_contigs):
        new_seq_with_most_overlap = min(new_contigs)
        # input_sequences = [seq for seq in input_sequences if seq not in new_seq_with_most_overlap]
        input_sequences.append(new_seq_with_most_overlap)

    print(max(input_sequences))
    return create_contig(input_sequences, original_seqs)


def combine_strings_with_overlap(sequence, check_sequence):
    for index, base in enumerate(sequence):
        sub_seq = sequence[index:]
        current_contig = ''
        for i, (sub_seq_base, check_seq_base) in enumerate(zip(sub_seq, check_sequence)):
            if sub_seq_base == check_seq_base:
                current_contig += sub_seq_base
                if i == len(sub_seq) - 1 and len(current_contig) > 1:
                    return f'{sequence[:-len(current_contig)]}{check_sequence}'
                    # First occurrence where overlap extends to the end of sub_seq will give
                    # combined string with most overlap. No need to keep checking.
            else:
                break


sequence_list = [x.get_seq() for x in create_fasta_list('rosalind_long.txt')]
genome = create_contig(sequence_list, sequence_list)
print(f'Longest superstring is {genome}')

with open('result.txt', 'w') as writer:
    writer.write(genome)
