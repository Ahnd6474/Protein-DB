from vae_module import Tokenizer, Config, load_vae, encode
from Bio import SeqIO
import re
cfg = Config(model_path="models/vae_epoch380.pt")
tok = Tokenizer.from_esm()

model = load_vae(cfg,
                 vocab_size=len(tok.vocab),
                 pad_idx=tok.pad_idx,
                 bos_idx=tok.bos_idx)

