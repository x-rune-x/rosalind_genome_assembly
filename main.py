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


class CombinedString:
    def __init__(self, sequence, overlap, parent1, parent2):
        self.sequence = sequence
        self.overlap = overlap
        self.parent1 = parent1
        self.parent2 = parent2


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


def create_superstring(input_sequences):
    if len(input_sequences) == 1:
        return input_sequences[0]

    new_string_combinations = []
    for sequence in input_sequences:
        # Want to find all other sequences with areas that overlap with the sequence we are checking.
        for check_sequence in input_sequences:
            if sequence != check_sequence:
                # Avoid checking duplicates.
                combo_string = combine_strings_with_overlap(sequence, check_sequence)
                if combo_string:
                    new_string_combinations.append(combo_string)

    new_string_combinations.sort(key=lambda x: x.overlap, reverse=True)
    new_seq_with_most_overlap = new_string_combinations[0]

    input_sequences.remove(new_seq_with_most_overlap.parent1)
    input_sequences.remove(new_seq_with_most_overlap.parent2)
    input_sequences.append(new_seq_with_most_overlap.sequence)

    return create_superstring(input_sequences)


def combine_strings_with_overlap(sequence, check_sequence):
    for index, base in enumerate(sequence):
        sub_seq = sequence[index:]
        current_overlap = ''
        for i, (sub_seq_base, check_seq_base) in enumerate(zip(sub_seq, check_sequence)):
            if sub_seq_base == check_seq_base:
                current_overlap += sub_seq_base
                if i == len(sub_seq) - 1 and len(current_overlap) > 1:
                    # First occurrence where overlap extends to the end of sub_seq will give
                    # combined string with most overlap. No need to keep checking.
                    combined_string = f'{sequence[:-len(current_overlap)]}{check_sequence}'
                    return CombinedString(combined_string, len(current_overlap), sequence, check_sequence)
            else:
                break


sequence_list = [x.get_seq() for x in create_fasta_list('rosalind_long.txt')]
genome = create_superstring(sequence_list)
print(f'Longest superstring is {genome}')

with open('result.txt', 'w') as writer:
    writer.write(genome)
