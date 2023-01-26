import Bio
from Bio import SeqIO
import sys
import re

fasta = sys.argv[1]
output = sys.argv[3]
target_chain = sys.argv[2]
chain_str = "Chain "+target_chain
auth_chain_str = "auth "+target_chain

handle = open(output,"w")
for record in SeqIO.parse(fasta,'fasta'):
    name = record.description
    # Here the case in which multiple chains are under the same fasta sequence
    multiple_chains = re.search('Chains([^|]*)\|.*', name, re.IGNORECASE)
    chains = multiple_chains.group(1) if multiple_chains else  ""
    if auth_chain_str in name or chain_str in name or target_chain in chains :
        SeqIO.write(record,output,"fasta")
