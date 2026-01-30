import numpy as np
<<<<<<< HEAD
import os
from tensorflow.keras.models import load_model
from qdrant_client import QdrantClient

# Remote Qdrant Cloud Client
client = QdrantClient(
    url="https://1603ecfd-8701-4b67-b638-e947cf4dc6cd.us-east4-0.gcp.cloud.qdrant.io", 
    api_key="eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJhY2Nlc3MiOiJtIn0.AAYBqPYFl9eIorZvc7vft3xLF7nXLVL4sPxijpNEXS8"
)

# Load model from the same directory as this script
MODEL_PATH = os.path.join(os.path.dirname(__file__), "Blog_simple_smi2lat.h5")

# Lazy loading to avoid crashing if model is missing
_model = None

def get_model():
    global _model
    if _model is None:
        print(f"🔍 Attempting to load SMILES model from: {MODEL_PATH}...")
        if os.path.exists(MODEL_PATH):
            try:
                _model = load_model(MODEL_PATH)
                print(f"✅ SUCCESS: SMILES model loaded successfully from {MODEL_PATH}")
            except Exception as e:
                print(f"❌ CRITICAL ERROR: Failed to load SMILES model: {e}")
        else:
            print(f"⚠️ WARNING: SMILES model file NOT FOUND at {MODEL_PATH}. Chemical search will be disabled.")
    return _model

# Trigger a load attempt when the module is first imported by main.py
# This ensures we see the status in the Docker logs immediately at startup.
get_model()

embed_dim = 28
charset = {'E', ']', '[', '(', '-', 'H', '4', '1', '#', 'O', ')', 'F', 'c', 'n', '+', 'o', 'C', '!', '=', 'N', '2', '3'}
char_to_int = {'E': 0, ']': 1, '[': 2, '(': 3, '-': 4, 'H': 5, '4': 6, '1': 7, '#': 8, 'O': 9, ')': 10, 'F': 11, 'c': 12, 'n': 13, '+': 14, 'o': 15, 'C': 16, '!': 17, '=': 18, 'N': 19, '2': 20, '3': 21}

def vectorize_single_smiles(smile, embed, char_to_int):
    one_hot = np.zeros((1, embed, len(char_to_int)), dtype=np.int8)
    one_hot[0, 0, char_to_int["!"]] = 1
    for j, c in enumerate(smile):
        if j + 1 >= embed: break # Truncate if too long
        if c in char_to_int:
            one_hot[0, j + 1, char_to_int[c]] = 1
    if len(smile) + 1 < embed:
        one_hot[0, len(smile) + 1:, char_to_int["E"]] = 1
    return one_hot[:, 0:-1, :]

def smiles_to_svg(smiles, size=(300, 300)):
    """Converts a SMILES string to a 2D SVG image string."""
    try:
        from rdkit import Chem
        from rdkit.Chem import Draw
        from rdkit.Chem.Draw import rdMolDraw2D
        
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
            
        d2d = rdMolDraw2D.MolDraw2DSVG(size[0], size[1])
        opts = d2d.drawOptions()
        opts.clearBackground = False # Better for dark themes
        
        d2d.DrawMolecule(mol)
        d2d.FinishDrawing()
        return d2d.GetDrawingText()
    except Exception as e:
        print(f"Drawing Error: {e}")
        return None

def search_similar_smiles(smiles_string, limit=5):
    model = get_model()
    if not model:
        return []
        
    try:
        X = vectorize_single_smiles(smiles_string, embed_dim, char_to_int)
        query_vec = model.predict(X)[0]
        
        # Use the remote cloud client and targeted collection
        results = client.query_points(
            collection_name="smiles_final",
            query=query_vec.tolist(),
            using="latent", # Ensure this matches your collection vector name
            limit=limit
        )
        
        return [
            {
                "smiles": r.payload.get("smiles", "Unknown"),
                "score": float(r.score),
                "svg": smiles_to_svg(r.payload.get("smiles", "")),
                "type": "molecule"
            }
            for r in results.points
        ]
    except Exception as e:
        print(f"[ERROR] SMILES Search failed: {e}")
        return []
=======
from tensorflow.keras.models import load_model
from qdrant_client import QdrantClient, models

client = QdrantClient(url = "https://1603ecfd-8701-4b67-b638-e947cf4dc6cd.us-east4-0.gcp.cloud.qdrant.io", api_key = "eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJhY2Nlc3MiOiJtIn0.AAYBqPYFl9eIorZvc7vft3xLF7nXLVL4sPxijpNEXS8")

smiles_to_latent_model = load_model("Blog_simple_smi2lat.h5")

embed = 28
charset = {'E', ']', '[', '(', '-', 'H', '4', '1', '#', 'O', ')', 'F', 'c', 'n', '+', 'o', 'C', '!', '=', 'N', '2', '3'}
char_to_int = {'E': 0, ']': 1, '[': 2, '(': 3, '-': 4, 'H': 5, '4': 6, '1': 7, '#': 8, 'O': 9, ')': 10, 'F': 11, 'c': 12, 'n': 13, '+': 14, 'o': 15, 'C': 16, '!': 17, '=': 18, 'N': 19, '2': 20, '3': 21}

def vectorize_single_smiles(smile, embed, charset, char_to_int):
    """
    Vectorize a single SMILES string for inference.
    
    Returns:
        X: shape (1, embed-1, len(charset))
    """
    one_hot = np.zeros((1, embed, len(charset)), dtype=np.int8)

    # Start token
    one_hot[0, 0, char_to_int["!"]] = 1

    # SMILES characters
    for j, c in enumerate(smile):
        if c not in char_to_int:
            raise ValueError(f"Unknown character '{c}' in SMILES")
        one_hot[0, j + 1, char_to_int[c]] = 1

    # End / padding token
    one_hot[0, len(smile) + 1:, char_to_int["E"]] = 1

    # Return only X (input to the model)
    return one_hot[:, 0:-1, :]

def search_similar_smiles(smiles_string):
    X = vectorize_single_smiles(smiles_string, embed, charset, char_to_int)
    query_vec = smiles_to_latent_model.predict(X)
    results = client.query_points(
        collection_name="smiles",
        query=query_vec,
        using="latent",
        limit=3
    )
    for point in results:
        print(point)
    return results
>>>>>>> d2ea12fb70eefdfe3e74f4e0253b9bb073cfc306
