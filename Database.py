from Protein import Protein

from Bio import SeqIO

def read_fasta(filename):
    data=[]
    for record in SeqIO.parse(filename, 'fasta'):
        data.append(Protein(record.id, record.description, record.locus, record.organism, record.seq))
    return data

class ProteinDB:
    def __init__(self,fasta=None,protein=None):
        if protein is None:
            protein=[]
        if fasta is None:
            fasta=[]
        else:
            read_fasta(fasta)
        self.proteins=fasta+protein
        self.proteins.sort()

    def __len__(self):
        return len(self.proteins)
    def __getitem__(self, index):
        return self.proteins[index]
    def __iter__(self):
        return iter(self.proteins)
    def __contains__(self, item):
        return item in self.proteins
    def __add__(self, other):
        self.proteins.extend(other.proteins)
        self.proteins=list(set(self.proteins))
        self.proteins.sort()
    def query(self,query):

