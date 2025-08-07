from protein_db import ProteinDB, Protein, ProteinQuery, generate_sequences


if __name__ == "__main__":
    # Vector-based search in simple mode
    db = ProteinDB(url="sqlite:///proteins_simple.db", simple=True)
    try:
        query1 = Protein("query", "", "", "", "MASTPFKFQLKGTINGKSFTVEGEGEGNSHEGSHKGKYVCTSGKLPMSWAALGTSFGYGMKYYTKYPSGLKNWFHEVMPEGFT")
        hits = db.search(query1.Z, k=5)
        print("Closest sequences:", hits)
        query2 = Protein("query", "", "", "", "MPRISDKLMKTRWRGFHSIPSIPPDLGGIYGIGEKTSRRKTTEHLYTGRAKDIKSRLMKHKYGHQAIDRKIRSNIKQKKLSDLRFKFVE")
    except Exception as exc:
        # Database may be empty if no sequences were loaded
        print(f"Vector search failed: {exc}")

    # Sequence generation using the genetic algorithm
    try:
        vectors = [query1.Z*(i/10)+query2.Z*(1-(i/10)) for i in range(1,10)]
        i=0
        for vector in vectors:
            i+=1
            seqs = generate_sequences(db, vector, cos_threshold=0.01, rmse_threshold=0.01)
            print(f"seq {i}:", seqs)
    except Exception as exc:
        print(f"Generation failed: {exc}")

'''    # Text-based queries on the full database
    db_full = ProteinDB(url="sqlite:///proteins.db")
    pq = ProteinQuery(db_full)
    for protein in pq.by_accession("ABC123"):
        print(protein.accession, protein.organism)
'''