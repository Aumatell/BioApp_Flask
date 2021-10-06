"""seqtools module
   contains functions for sequence cleaning and visualization"""
    

def clean(seq):
    """Returns a sequence (seq) without any space or end-of-line character"""
    return seq.strip().replace(' ', '').replace('\n', '').replace('\r', '')


def chunks(seq, size):
    """Returns a sequence (seq) formatted in chuncks of size (size) sepparated
       with end-of-line characters"""
    seq = str(seq)  # in case seq is a Biopython sequence
    seq_list = [seq[i:i+size] for i in range(0, len(seq), size)]
    return '\n'.join(seq_list)


def cutseq(seq):
    """Returns a sequence that is multiple of 3 by cutting the last bases.
       This is useful in translation"""
    rem = len(seq) % 3
    if rem != 0:
        return seq #[:-rem]
    else:
        return seq


# This will only be executed when the file is executed as a script
# Not when it is imported
if __name__ == '__main__':
    with open('test_sequence1.txt') as file_input:
        seq = file_input.read()
    print(chunks(clean(seq), 100))

    with open('test_sequence2.txt') as file_input:
        seq = file_input.read()
    print(chunks(clean(seq), 100))

  
