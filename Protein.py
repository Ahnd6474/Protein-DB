class Protein:
    def __init__(self,seq_id,seq,lat):
        self.seq_id = seq_id
        self.sequence = seq
        self.latent_vec=lat
    def __len__(self):
        return len(self.sequence)
    def __getitem__(self,idx):
        return self.sequence[idx]

class 

def comb_sort(l, s):
    n = len(l)
    gap = n
    while True:
        gap = int(gap // s) if gap // s > 1 else 1
        c = 0 if gap == 1 else 1
        for i in range(n - gap):
            if l[i] > l[i + gap]:
                l[i], l[i + gap] = l[i + gap], l[i]
                c += 1
        if c == 0:
            break
    return l

def insertion_sort(l):
    n=len(l)
    for i in range(n):
        t=l[i]
        c=i
        while c>0 and l[c-1]>=t:
            l[c]=l[c-1]
            c-=1
        l[c]=t
    return l
