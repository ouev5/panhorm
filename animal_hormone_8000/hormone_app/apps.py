# hormone_app/apps.py
from django.apps import AppConfig
import logging

# Configure logger
logger = logging.getLogger(__name__)

class HormoneAppConfig(AppConfig):
    default_auto_field = 'django.db.models.BigAutoField'
    name = 'hormone_app'

    def ready(self):
        """Initialize vector data when Django app starts"""
        try:
            # Lazy import to avoid circular dependencies
            from .utils.rag_retriever import vector_store
            # Load database data into vector index
            vector_store.load_data_to_vector()
            logger.info("? Hormone data vector index loaded successfully")
        except Exception as e:
            logger.error(f"? Failed to load vector index: {str(e)}", exc_info=True)
