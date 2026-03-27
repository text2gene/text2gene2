"""CLI: hgvs2pmid"""
import asyncio
import json
import sys

import typer

app = typer.Typer(help="text2gene2 — find PubMed literature for an HGVS variant")


@app.command()
def search(
    hgvs: str = typer.Argument(..., help="HGVS variant string"),
    output: str = typer.Option("table", "--output", "-o", help="Output format: table | json"),
):
    """Search all sources for a single variant."""
    from text2gene2.pipeline import query_variant

    result = asyncio.run(query_variant(hgvs))

    if output == "json":
        typer.echo(result.model_dump_json(indent=2))
        return

    typer.echo(f"\nVariant:  {result.input_hgvs}")
    if result.lvg and result.lvg.gene_symbol:
        typer.echo(f"Gene:     {result.lvg.gene_symbol}")
    if result.lvg and result.lvg.rsids:
        typer.echo(f"rsIDs:    {', '.join(result.lvg.rsids)}")

    typer.echo(f"\nFound {len(result.citations)} unique PMIDs:\n")
    for citation in result.citations:
        src_str = " ".join(f"[{s}]" for s in citation.sources)
        typer.echo(f"  PMID {citation.pmid:>10}  conf={citation.confidence}  {src_str}")

    if result.errors:
        typer.echo("\nSource errors:")
        for src, err in result.errors.items():
            typer.echo(f"  {src}: {err}", err=True)


@app.command()
def batch(
    file: typer.FileText = typer.Argument(..., help="File with one HGVS per line"),
    output: str = typer.Option("json", "--output", "-o", help="Output format: json | tsv"),
):
    """Search all sources for a file of variants (one HGVS per line)."""
    from text2gene2.pipeline import query_variant

    variants = [line.strip() for line in file if line.strip() and not line.startswith("#")]
    typer.echo(f"Processing {len(variants)} variants...", err=True)

    async def run_all():
        tasks = [query_variant(v) for v in variants]
        return await asyncio.gather(*tasks, return_exceptions=True)

    results = asyncio.run(run_all())

    if output == "json":
        out = []
        for v, r in zip(variants, results):
            if isinstance(r, Exception):
                out.append({"variant": v, "error": str(r)})
            else:
                out.append(r.model_dump())
        typer.echo(json.dumps(out, indent=2))
    else:
        typer.echo("variant\tpmid\tsources\tconfidence")
        for v, r in zip(variants, results):
            if isinstance(r, Exception):
                typer.echo(f"{v}\tERROR\t{r}\t0")
            else:
                for c in r.citations:
                    typer.echo(f"{v}\t{c.pmid}\t{'|'.join(c.sources)}\t{c.confidence}")


@app.command()
def serve(
    host: str = typer.Option("0.0.0.0", help="Bind host"),
    port: int = typer.Option(8000, help="Bind port"),
    reload: bool = typer.Option(False, "--reload", help="Auto-reload (dev only)"),
):
    """Start the web server."""
    import uvicorn
    uvicorn.run("text2gene2.api.app:app", host=host, port=port, reload=reload)


if __name__ == "__main__":
    app()
