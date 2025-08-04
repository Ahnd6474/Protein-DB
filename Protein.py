from vae_module import Tokenizer, Config, load_vae, encode
from Bio import SeqIO
import re
cfg = Config(model_path="models/vae_epoch380.pt")
tok = Tokenizer.from_esm()

model = load_vae(cfg,
                 vocab_size=len(tok.vocab),
                 pad_idx=tok.pad_idx,
                 bos_idx=tok.bos_idx)

class Protein:
    def __init__(self,accession, description, locus, organism,seq):
        self.accession = accession
        self.description = description
        self.locus = locus
        self.organism = organism
        self.sequence = seq
        self.Z=np.array(encode(model, seq, tok, cfg.max_len))
    def __len__(self):
        return len(self.sequence)
    def __getitem__(self,idx):
        return self.sequence[idx]
    def __sub__(self,other):
        return abs(self.Z-other.Z)
    def __repr__(self):
        return str({
            "accession": self.accession,
            "description": self.description,
            "locus": self.locus,
            "organism": self.organism,
            "sequence": self.sequence,
        })
    def __str__(self):
        return str({
            "accession": self.accession,
            "description": self.description,
            "locus": self.locus,
            "organism": self.organism,
            "sequence": self.sequence,
        })

def describe(description):
    # accession (첫 번째 필드)
    parts = description.split(None, 1)
    accession = parts[0]
    rest = parts[1] if len(parts) > 1 else ""

    # 유기체명 (대괄호 안)
    org_match = re.search(r"\[([^\]]+)\]$", rest)
    organism = org_match.group(1) if org_match else ""

    # 유기체명 제거
    rest = re.sub(r"\s*\[[^\]]+\]$", "", rest)

    # 나머지 설명과 로커스 분리
    # 예시: "hypothetical protein TI39_contig5958g00003"
    desc_parts = rest.rsplit(None, 1)
    if len(desc_parts) == 2:
        description, locus = desc_parts
    else:
        description, locus = rest, ""

    return accession,description,locus,organism

