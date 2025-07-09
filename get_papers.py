# get_papers.py
import argparse
import sys
from pubmed_fetcher.fetcher import PubMedFetcher

def main():
    parser = argparse.ArgumentParser(
        description="Fetch research papers from PubMed with non-academic authors from pharma/biotech companies"
    )
    parser.add_argument("query", help="PubMed query string")
    parser.add_argument("-d", "--debug", action="store_true", help="Enable debug logging")
    parser.add_argument("-f", "--file", help="Output file (CSV format)")
    
    args = parser.parse_args()
    
    # Initialize fetcher
    fetcher = PubMedFetcher(
        email="shreyasinghin24@gmail.com",  # add your email id to access PubMed API
        debug=args.debug
    )
    
    # Fetch papers
    papers = fetcher.fetch_papers(args.query)
    
    # Write results
    fetcher.write_csv(papers, args.file)

if __name__ == "__main__":
    main()