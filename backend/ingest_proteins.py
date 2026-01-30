import requests
import json
import os
from qdrant_db import init_collection, PROTEIN_COLLECTION, client, PointStruct
from embeddings import embed
from datetime import datetime

def fetch_protein_library(limit=50):
    """Scrapes protein metadata from RCSB PDB."""
    search_url = "https://search.rcsb.org/rcsbsearch/v2/query"
    
    query = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_entity_source_organism.scientific_name",
                "operator": "exact_match",
                "value": "Homo sapiens"
            }
        },
        "request_options": {
            "paginate": {"start": 0, "rows": limit},
            "sort": [{"sort_by": "rcsb_accession_info.deposit_date", "direction": "desc"}]
        },
        "return_type": "entry"
    }

    print(f"📡 Requesting {limit} Human structures from PDB...")
    
    try:
        response = requests.post(search_url, json=query)
        response.raise_for_status()
        results = response.json().get("result_set", [])
        pdb_ids = [res["identifier"] for res in results]
    except Exception as e:
        print(f"❌ Search failed: {e}")
        return []

    final_library = []

    for pdb_id in pdb_ids:
        try:
            print(f"🧬 Processing: {pdb_id}...", end="\r")
            data_url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
            res = requests.get(data_url)
            if res.status_code != 200: continue
            
            data = res.json()
            raw_title = data.get("struct", {}).get("title", "Unknown Protein")
            clean_name = raw_title.strip().capitalize()
            
            desc = data.get("struct_keywords", {}).get("pdbx_details")
            if not desc:
                desc = data.get("struct_keywords", {}).get("text", "Molecular Structure")

            stoich = data.get("rcsb_entry_info", {}).get("polymer_composition", "Unknown")

            final_library.append({
                "pdb_id": pdb_id,
                "name": clean_name,
                "organism": "Homo sapiens",
                "description": desc,
                "stoichiometry": stoich
            })
        except Exception as e:
            print(f"\n⚠️ Error on {pdb_id}: {e}")

    return final_library

def run_protein_ingestion():
    print("[INFO] Starting Protein Ingestion Pipeline...")
    init_collection(PROTEIN_COLLECTION)
    
    proteins = fetch_protein_library(limit=30) # Start with 30 for testing
    if not proteins:
        print("[ERROR] No proteins fetched.")
        return

    print(f"\n✅ Fetched {len(proteins)} proteins. Generating BGE embeddings...")
    
    points = []
    for idx, protein in enumerate(proteins):
        pdb_id = protein["pdb_id"]
        # Creating a rich searchable string as requested
        searchable_text = f"PDB ID: {pdb_id}. Name: {protein['name']}. Details: {protein['description']}. Stoichiometry: {protein['stoichiometry']}"
        
        # Using unified BGE Small embedder
        vector = embed([searchable_text])[0].tolist()
        
        points.append(PointStruct(
            id=idx, # Simple numerical ID for this collection
            vector={"text": vector}, # Must match VECTOR_NAME in qdrant_db
            payload=protein
        ))

    print(f"🚀 Uploading {len(points)} proteins to '{PROTEIN_COLLECTION}'...")
    client.upsert(collection_name=PROTEIN_COLLECTION, points=points)
    print("✨ Protein vector engine ready!")

if __name__ == "__main__":
    run_protein_ingestion()
