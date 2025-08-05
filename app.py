"""Interactive Streamlit app for querying the protein database."""

import streamlit as st

from protein_db.database import ProteinDB
from protein_db.query import ProteinQuery
from protein_db.blast import deep_blast
from protein_db.visualize import plot_embeddings


st.title("Protein Database UI")

sequence = st.text_area("Query sequence")
threshold = st.number_input("Threshold", value=0.5)
metric = st.selectbox("Metric", ["euclidean", "cosine"])
db_url = st.text_input("Database URL", "sqlite:///proteins.db")

if st.button("Search") and sequence:
    db = ProteinDB(db_url)
    query = ProteinQuery(db)
    hits = query.by_embedding(sequence, threshold, metric)
    st.write(f"{len(hits)} hits found")

    if hits:
        fig = plot_embeddings(hits, sequence)
        st.pyplot(fig)

        if st.checkbox("Run BLAST on hits"):
            results = deep_blast(sequence, db, threshold, metric)
            for accession, result in results:
                st.subheader(accession)
                st.text(result)
