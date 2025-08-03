from Protein import Protein

from Bio import SeqIO

def read_fasta(filename):
    data=[]
    for record in SeqIO.parse(filename, 'fasta'):
        data.append(Protein(record.id, record.description, record.locus, record.organism, record.seq))
    return data
class ProteinDB:
    def __init__(self):

