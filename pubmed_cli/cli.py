import click 
from pubmed_cli.pubmed import fetch_pubmed
from pubmed_cli.filters import is_non_academic
from pubmed_cli.utils import save_papers


@click.group()
def cli():
    """PubMed CLI Tool"""
    pass

@cli.command()
@click.argument("query")
@click.option("--max", "max_results", default=100, help="Number of papers to be fetched")
@click.option("--file", "output_file", type=click.Path(), help="Export results to SV file")
@click.option("--debug", is_flag=True, help="Enable debug mode")
def search(query, max_results, output_file, debug):
    """ Search PubMed with a query and advanced options/filters."""
    if debug:
        click.echo(f" Searching PubMed for: '{query}' (max={max_results})")

    papers = fetch_pubmed(query, max_results=max_results, debug=debug)
    if debug:
        click.echo(f" Found {len(papers)} papers")

    filtered_papers = []
    for paper in papers:
        affiliations = paper.get("Affiliations", [])
        if any(is_non_academic(aff) for aff in affiliations):
            filtered_papers.append(paper)
  
    if output_file:
        save_papers(filtered_papers, output_file)
        click.echo(f"Results saved, Exported{len(filtered_papers)} papers to {output_file}")
    else:
        for paper in filtered_papers:
            click.echo("---------")
            click.echo(f" Title: {paper['Title']}")
            click.echo(f" Affiliations: {', '.join(paper['Affiliations'])}")
    

def hello():
    """Welcome message for the CLI"""
    click.echo("Hello from PubMed CLI!")

if __name__ == "__main__":
    cli()
