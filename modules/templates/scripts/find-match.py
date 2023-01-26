import sys
from Bio import SeqIO

target = sys.argv[1]
query = sys.argv[2]


def get_seq(input_file):
    fasta_sequences = SeqIO.parse(open(input_file),'fasta')
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        return(sequence)

def posit_in_seq(string,site):

    len_site = len(site)
    if site in string:
        for i in range (0, len(string)):
            if site == string[i:i+len_site]:
                return i+1
    elif site not in string:
        return -1

print(posit_in_seq(get_seq(target), get_seq(query)))
