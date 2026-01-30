import numpy as np
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
