"""Versioned entry point for the manuscript GeneReg v2.6 application."""
from app import app, _do_analysis, HTML_PAGE, gf, ca, ma, la, ai_a

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8003)
