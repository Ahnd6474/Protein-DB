from vae_module import Tokenizer, Config, load_vae, encode

cfg = Config(model_path="models/vae_epoch380.pt")
tok = Tokenizer.from_esm()

model = load_vae(cfg,
                 vocab_size=len(tok.vocab),
                 pad_idx=tok.pad_idx,
                 bos_idx=tok.bos_idx)

class Protein:
    def __init__(self,seq_id,seq):
        self.seq_id = seq_id
        self.sequence = seq
        self.Z=encode(model, seq, tok, cfg.max_len)
    def __len__(self):
        return len(self.sequence)
    def __getitem__(self,idx):
        return self.sequence[idx]

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
