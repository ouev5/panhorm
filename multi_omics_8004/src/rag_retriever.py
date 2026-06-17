import os
import glob
from typing import List, Dict, Any, Optional
from langchain_core.documents import Document


class RAGRetriever:
    def __init__(self, persist_directory: str = "./rag_index"):
        self.persist_directory = persist_directory
        self.vectorstore = None
        self.bm25_retriever = None
        self.embeddings = None  # Lazy loaded
        self.text_splitter = None

    def _ensure_deps(self):
        """Lazy load heavy dependencies only when needed."""
        if self.text_splitter is None:
            from langchain_text_splitters import RecursiveCharacterTextSplitter
            self.text_splitter = RecursiveCharacterTextSplitter(
                chunk_size=1000, chunk_overlap=200, length_function=len, is_separator_regex=False,
            )
        if self.embeddings is None:
            from langchain_community.embeddings import HuggingFaceEmbeddings
            self.embeddings = HuggingFaceEmbeddings(
                model_name="BAAI/bge-large-zh",
                model_kwargs={"device": "cpu"},
                encode_kwargs={"normalize_embeddings": True}
            )

    def _load_document(self, file_path: str) -> List[Document]:
        try:
            if file_path.endswith(".pdf"):
                from langchain_community.document_loaders import PyPDFLoader
                loader = PyPDFLoader(file_path)
            elif file_path.endswith(".md"):
                from langchain_community.document_loaders import UnstructuredMarkdownLoader
                loader = UnstructuredMarkdownLoader(file_path)
            else:
                from langchain_community.document_loaders import TextLoader
                loader = TextLoader(file_path, encoding="utf-8")
            return loader.load()
        except Exception as e:
            print(f"Error loading document {file_path}: {e}")
            return []

    def _load_documents_from_directory(self, directory_path: str) -> List[Document]:
        documents = []
        for ext in [".pdf", ".md", ".txt", ".text"]:
            for fpath in glob.glob(os.path.join(directory_path, f"**/*{ext}"), recursive=True):
                print(f"Loading document: {fpath}")
                documents.extend(self._load_document(fpath))
        return documents

    def build_index(self, documents_path: str) -> Dict[str, Any]:
        self._ensure_deps()
        if os.path.isdir(documents_path):
            documents = self._load_documents_from_directory(documents_path)
        else:
            documents = self._load_document(documents_path)
        chunks = self.text_splitter.split_documents(documents)
        from langchain_community.vectorstores import Chroma
        self.vectorstore = Chroma.from_documents(
            documents=chunks, embedding=self.embeddings, persist_directory=self.persist_directory
        )
        from langchain_community.retrievers import BM25Retriever
        self.bm25_retriever = BM25Retriever.from_documents(chunks)
        if self.vectorstore:
            self.vectorstore.persist()
        return {
            "total_documents": len(documents), "total_chunks": len(chunks),
            "persist_directory": self.persist_directory, "embedding_model": "BAAI/bge-large-zh",
            "chunk_size": 1000, "chunk_overlap": 200
        }

    def retrieve(self, query: str, top_k: int = 5) -> List[Document]:
        self._ensure_deps()
        if os.path.exists(self.persist_directory):
            from langchain_community.vectorstores import Chroma
            self.vectorstore = Chroma(
                persist_directory=self.persist_directory, embedding_function=self.embeddings
            )
            return self.vectorstore.similarity_search(query, k=top_k)
        raise ValueError("No index built. Please call build_index first.")

    def load_index(self, persist_directory: Optional[str] = None):
        self._ensure_deps()
        if persist_directory:
            self.persist_directory = persist_directory
        if os.path.exists(self.persist_directory):
            from langchain_community.vectorstores import Chroma
            self.vectorstore = Chroma(
                persist_directory=self.persist_directory, embedding_function=self.embeddings
            )
            print(f"Loaded index from {self.persist_directory}")
        else:
            raise ValueError(f"Index directory not found: {self.persist_directory}")
