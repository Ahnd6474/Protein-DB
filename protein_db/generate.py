import random
from typing import List, Sequence

import numpy as np

from .database import ProteinDB
from .protein import Protein

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"

def cosine_similarity(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=np.float32)
    b = np.asarray(b, dtype=np.float32)
    return float(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b)))

def rmse(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=np.float32)
    b = np.asarray(b, dtype=np.float32)
    return float(np.sqrt(np.mean((a - b) ** 2)))

def recombine(a: str, b: str) -> str:
    cut = random.randint(1, min(len(a), len(b)) - 1)
    return a[:cut] + b[cut:]

def substitute(seq: str) -> str:
    if not seq:
        return seq
    idx = random.randrange(len(seq))
    aa = random.choice(AMINO_ACIDS)
    return seq[:idx] + aa + seq[idx + 1 :]

def insert(seq: str) -> str:
    idx = random.randrange(len(seq) + 1)
    aa = random.choice(AMINO_ACIDS)
    return seq[:idx] + aa + seq[idx:]

def delete(seq: str) -> str:
    if not seq:
        return seq
    idx = random.randrange(len(seq))
    return seq[:idx] + seq[idx + 1 :]

def mutate(seq: str) -> str:
    op = random.random()
    if op < 0.25:
        return substitute(seq)
    elif op < 0.5:
        return insert(seq)
    elif op < 0.75:
        return delete(seq)
    else:
        return seq  # no-op fallback

def generate_sequences(
    db: ProteinDB,
    vector: Sequence[float],
    cos_threshold: float = 0.1,
    rmse_threshold: float = 0.1,
    generations: int = 10,
) -> List[str]:
    """Generate sequences close to ``vector`` using a simple genetic algorithm.

    Parameters
    ----------
    db:
        ``ProteinDB`` instance in simple mode.
    vector:
        Target 256-dimensional embedding vector.
    cos_threshold:
        Maximum allowed cosine distance for accepting mutations.
    rmse_threshold:
        Maximum allowed RMSE for accepting mutations.
    generations:
        Number of generations to evolve.

    Returns
    -------
    List[str]
        Five sequences most similar to ``vector``.
    """

    target = np.asarray(vector, dtype=np.float32)
    population = db.search(target, k=100)
    if len(population) < 100:
        population = list(population)
    else:
        population = list(population)
    while len(population) < 200:
        population.append(mutate(random.choice(population)))

    for _ in range(generations):
        scored = []
        for seq in population:
            emb = Protein("tmp", "", "", "", seq).Z
            cos = cosine_similarity(target, emb)
            error = rmse(target, emb)
            scored.append((cos, error, seq))
        scored.sort(key=lambda x: (-x[0], x[1]))
        elites = [s[2] for s in scored[:50]]
        next_pop = elites.copy()
        while len(next_pop) < 200:
            op = random.random()
            if op < 0.25 and len(population) >= 2:
                s1, s2 = random.sample(population, 2)
                child = recombine(s1, s2)
            else:
                child = mutate(random.choice(population))
            emb = Protein("tmp", "", "", "", child).Z
            cos = cosine_similarity(target, emb)
            error = rmse(target, emb)
            if (1 - cos) <= cos_threshold and error <= rmse_threshold:
                next_pop.append(child)
        population = next_pop

    final_scores = []
    for seq in population:
        emb = Protein("tmp", "", "", "", seq).Z
        cos = cosine_similarity(target, emb)
        error = rmse(target, emb)
        final_scores.append((cos, error, seq))
    final_scores.sort(key=lambda x: (-x[0], x[1]))
    return [s[2] for s in final_scores[:5]]
