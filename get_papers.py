import argparse
from pubmed_fetcher.fetcher import PubMedFetcher

def get_papers_list():
    """Entry point for poetry script"""
    parser = argparse.ArgumentParser(
        description="Fetch research papers from PubMed with non-academic authors from pharma/biotech companies"
    )
    parser.add_argument("query", help="PubMed query string")
    parser.add_argument("-d", "--debug", action="store_true", help="Enable debug logging")
    parser.add_argument("-f", "--file", help="Output file (CSV format)")
    
    args = parser.parse_args()
    
    fetcher = PubMedFetcher(
        email="shreyasinghin24@gmail.com",
        debug=args.debug
    )
    
    papers = fetcher.fetch_papers(args.query)
    fetcher.write_csv(papers, args.file)

if __name__ == "__main__":
    get_papers_list()