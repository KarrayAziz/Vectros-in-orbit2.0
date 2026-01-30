from fastapi import FastAPI
from ingest import update_database
from smiles import search_similar_smiles
from qdrant_client import QdrantClient
import os
from tensorflow.keras.models import load_model
from qdrant_client import QdrantClient, models

from fastapi.middleware.cors import CORSMiddleware

app = FastAPI()

# Add CORS middleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # In production, specify your frontend URL
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Read Qdrant host/port from environment so Docker Compose service name works
QDRANT_HOST = os.getenv("QDRANT_HOST", "localhost")
QDRANT_PORT = int(os.getenv("QDRANT_PORT", 6333))

client = QdrantClient(host=QDRANT_HOST, port=QDRANT_PORT)


@app.get("/")
def root():
    return {"status": "ok", "message": "Bio-Vector backend running"}


@app.post("/update-db")
def update_db():
    return update_database()


@app.get("/search")
def semantic_search(query: str, limit: int = 12, offset: int = 0, search_type: str = "text"):
    from embeddings import embed
    from qdrant_db import COLLECTION, PROTEIN_COLLECTION
    from smiles import search_similar_smiles

    
    # Determine collection name based on search type
    if search_type == "molecule":
        return search_similar_smiles(query, limit=limit)

    vector = embed([query])[0]
    
    collection_name = PROTEIN_COLLECTION if search_type == "protein" else COLLECTION

    response = client.query_points(
        collection_name=collection_name,
        query=vector.tolist(),
        using="text",
        limit=limit,
        offset=offset
    )
    results = response.points

    return [
        {
            "pmid": r.payload["pmid"],
            "title": r.payload["title"],
            "abstract": r.payload.get("chunk_text", r.payload.get("abstract", "")),  # Return chunk_text if available
            "url": r.payload["url"],
            "score": r.score
        }
        for r in results
    ]

@app.get("/smiles-search")
def smiles_search(smiles: str, limit: int = 3):
    results = search_similar_smiles(smiles)
    return results
