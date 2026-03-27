from fastapi import FastAPI
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
from pathlib import Path

from .routes import router


def create_app() -> FastAPI:
    app = FastAPI(
        title="text2gene2",
        description="Given an HGVS genetic variant, find all relevant PubMed literature.",
        version="2.0.0",
        docs_url="/api/docs",
        redoc_url="/api/redoc",
    )

    app.include_router(router)

    # Static files (CSS, JS)
    static_dir = Path(__file__).parent.parent / "static"
    if static_dir.exists():
        app.mount("/static", StaticFiles(directory=static_dir), name="static")

    return app


app = create_app()
