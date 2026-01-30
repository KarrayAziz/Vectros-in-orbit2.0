# Bio-Vector Orbit - Discovery Engine v2.0
<img width="1918" height="902" alt="image" src="https://github.com/user-attachments/assets/9ba6a76b-56fc-4502-beb0-92cfd337965c" />

A **Vector-Powered Discovery Engine** that transforms fragmented biological data (papers, DNA sequences, and chemical molecules) into a unified, searchable intelligence layer using **NCBI**, **Qdrant**, **FastEmbed**, and **React**.

![Bio-Vector Orbit](https://img.shields.io/badge/Status-Hackathon%20Ready-success)
![Docker](https://img.shields.io/badge/Docker-Enabled-blue)
![License](https://img.shields.io/badge/License-MIT-green)

---

## 🚀 Quick Start with Docker (Recommended for Judges)

**Prerequisites:** [Docker](https://www.docker.com/get-started) and Docker Compose installed.

### 1. Start the Application

```bash
docker-compose up
```

That's it! The entire stack will start automatically.

### 2. Access the Application

- **Frontend Dashboard**: [http://localhost:3000](http://localhost:3000)
- **Backend API**: [http://localhost:8000/docs](http://localhost:8000/docs) (Interactive Swagger UI)
- **Qdrant Dashboard**: [http://localhost:6333/dashboard](http://localhost:6333/dashboard)

### 3. Try It Out

1. Open [http://localhost:3000](http://localhost:3000)
2. Type "**Insulin**" or "**BRCA1**" in the search bar
3. Hit **Enter**
4. Watch as the system:
   - Fetches real papers from NCBI PubMed
   - Processes them with semantic chunking
   - Displays scientifically relevant results
   - Provides links back to original sources

### 4. Stop the Application

```bash
docker-compose down
```

---

## ✨ Features

### Core Capabilities
- 🧬 **Real-Time NCBI Search**: Automatically fetches the latest papers from PubMed on each query
- 🔍 **Semantic Vector Search**: Uses FastEmbed (BAAI/bge-small-en-v1.5) for meaning-based retrieval
- 📊 **Physics-Informed Ranking**: Results include thermodynamic stability (ΔG values)
- 🧱 **Semantic Chunking**: Chonkie library breaks abstracts at logical boundaries
- 🏗️ **Multimodal Embeddings**: Supports Text, Protein, and Molecule vectors
- 🔬 **3D Visualization**: iCn3D integration for protein structure viewing
- 📎 **Scientific Traceability**: Every result links back to the original PubMed article

### Technical Highlights
- **Pagination Support**: Seamlessly load more results (12 per request) to explore deep into the search space.
- **Named Vectors**: Separate vector spaces for Text (384D), Protein (384D), Molecule (384D)
- **HNSW Indexing**: High-precision clustering (m=16, ef_construct=200) for fast retrieval
- **MMR Diversity**: Maximal Marginal Relevance for diverse result sets
- **Zero-Setup**: Full Docker orchestration handles all dependencies and initialization automatically.

---

## 🏗️ Architecture

```
┌─────────────┐      ┌──────────────┐      ┌─────────────┐
│   React     │─────▶│   FastAPI    │─────▶│   Qdrant    │
│  Frontend   │      │   Backend    │      │  Vector DB  │
│  (Port 80)  │      │  (Port 8000) │      │ (Port 6333) │
└─────────────┘      └──────────────┘      └─────────────┘
                            │
                            ▼
                     ┌──────────────┐
                     │ NCBI PubMed  │
                     │   (Biopython)│
                     └──────────────┘
```

### Technology Stack
- **Frontend**: React 19 + Vite + TypeScript + TailwindCSS
- **Backend**: Python 3.1 + FastAPI + Uvicorn
- **Vector Database**: Qdrant (with Named Vectors)
- **Embeddings**: FastEmbed (BAAI/bge-small-en-v1.5)
- **Chunking**: Chonkie (Semantic Chunker)
- **Data Source**: NCBI PubMed (via Biopython Entrez)
- **3D Visualization**: iCn3D

---

## 📁 Project Structure

```
Vectros-in-orbit2.0/
├── backend/
│   ├── main.py              # FastAPI application & API endpoints
│   ├── embeddings.py        # FastEmbed model management
│   ├── pubmed.py            # NCBI PubMed retrieval logic
│   ├── qdrant_db.py         # Qdrant client & vector operations
│   ├── ingest.py            # Data pipeline & semantic chunking
│   ├── create_collection.py # Collection initialization script
│   ├── requirements.txt     # Python dependencies
│   └── Dockerfile           # Backend container definition
├── components/
│   ├── MoleculeViewer.tsx   # iCn3D 3D molecular viewer
│   ├── ResultCard.tsx       # Interactive search result card
│   └── ...                  # Other UI components
├── App.tsx                  # Main React application & state management
├── constants.ts             # Application constants & configuration
├── docker-compose.yml       # Multi-container orchestration
├── Dockerfile               # Production frontend container (Nginx)
└── README.md                # This documentation
```

---

## 🧪 How It Works

1. **User Query**: User types a scientific query (e.g., "Insulin")
2. **NCBI Fetch**: Backend queries PubMed via Biopython to grab the most recent relevant abstracts.
3. **Semantic Chunking**: Chonkie splits abstracts intelligently at semantic boundaries.
4. **Vector Embedding**: FastEmbed converts each chunk into a high-dimensional vector.
5. **Qdrant Storage**: Vectors + rich metadata (PMID, title, source, chunk text) are upserted into Qdrant.
6. **Similarity Search**: The query is embedded and compared against the vector store using Cosine similarity.
7. **Pagination**: Users can click "Display next 12 results" to increment the search offset and load deeper insights.
8. **Visualization**: One-click 3D structure viewing for proteins via iCn3D.

---

---

## 🐳 Running with Docker

The application is fully containerized for a seamless experience across all platforms.

```bash
# Build and start all services
docker-compose up --build

# Shutdown the application
docker-compose down
```

---

## 🎯 Use Cases

- **Drug Discovery**: Find papers on specific protein targets or compounds
- **Literature Review**: Semantic search across biomedical abstracts
- **Hypothesis Generation**: Discover connections between biological concepts
- **Education**: Explore biological topics with visual 3D structures

---

## 📝 License

MIT License - See LICENSE file for details

---

## 👥 Author

Built for the **Vectors in Orbit** Hackathon by **Karray Aziz**

---

## 🙏 Acknowledgments

- **NCBI PubMed**: For providing open access to biomedical literature
- **Qdrant**: For the powerful vector database solution
- **FastEmbed**: For efficient, lightweight embeddings
- **iCn3D**: For 3D molecular visualization

---

**Questions?** Open an issue or contact the maintainer.

**Enjoy exploring the biology discovery engine!** 🧬🚀
