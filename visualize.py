"""Visualization utilities for protein embeddings."""

from typing import List

import matplotlib.pyplot as plt
import numpy as np
from sklearn.decomposition import PCA

from Database import ProteinModel
from Protein import Protein


def plot_embeddings(proteins: List[ProteinModel], query_sequence: str) -> plt.Figure:
    """Plot protein embeddings with the query highlighted.

    Parameters
    ----------
    proteins:
        List of ``ProteinModel`` objects to include in the plot.
    query_sequence:
        Amino acid sequence used for the query.

    Returns
    -------
    matplotlib.figure.Figure
        Figure object containing a 2D PCA scatter plot of embeddings.
    """

    query_embedding = Protein("query", "", "", "", query_sequence).Z
    embeddings = [p.embedding for p in proteins] + [query_embedding]
    pca = PCA(n_components=2)
    reduced = pca.fit_transform(np.vstack(embeddings))

    fig, ax = plt.subplots()
    ax.scatter(reduced[:-1, 0], reduced[:-1, 1], label="proteins")
    ax.scatter(reduced[-1, 0], reduced[-1, 1], color="red", label="query")
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title("Protein Embeddings")
    ax.legend()
    return fig
