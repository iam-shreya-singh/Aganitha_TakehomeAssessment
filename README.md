
#  PubMed Pharma/Biotech Paper Fetcher

A Python program to fetch research papers from PubMed with non-academic authors from pharmaceutical or biotech companies.


## Features

-Fetch papers using PubMed's full query syntax

-Filter for papers with ≥1 author from pharma/biotech companies

-Output results as CSV with required columns

-Command-line interface with debug mode and file output options

## Installation

Install Poetry: https://python-poetry.org/docs/

Install dependencies:

```bash
 poetry install
```

## Dependencies
Python > 3.8+

BioPython (for PubMed API access)

Poetry (dependency management)


## Project Structure

```bash
project/
├── pubmed_fetcher/    # ← MODULE (backend logics)
│   ├── __init__.py
│   └── fetcher.py     # ← Contains PubMedFetcher class
├── get_papers.py      # ← CLI PROGRAM (user interface)
├── pyproject.toml
└── README.md
```


## Usage/Examples

```bash
poetry run get-papers-list "your query"

poetry run get-papers-list "cancer AND (vaccine OR immunotherapy) AND 2020-2023[PDAT]"

poetry run get-papers-list "cancer AND (vaccine OR immunotherapy) AND company[AFFILIATION] AND 2020-2023[PDAT]"

```


## Options

-d, --debug: Enable debug logging

-f, --file: Save your results to CSV file

-h, --help: Show help message

### Save to file

```bash
poetry run get-papers-list "query" --file results.csv
```


## How It Works

 1. Uses BioPython's Entrez module to access PubMed API
 2. Parses author affiliations to identify non-academic institutions
 3. Filters for pharma/biotech companies usingkeyword matching
 4. Extracts corresponding author email (first author with "corresponding" in affiliation)
 5. Outputs results in specified CSV format


### Heuristics for Non-academic Authors

--Pharma/Biotech: "pharma", "pharmaceutical", "biotech", "biotechnology", "inc", "llc", "ltd", "sa", "ag", "co", "corp", " Holding", "Holdings", " GmbH", " AB", " S.A.", "Pharmaceuticals", "Biotechnology", "Bioscience", "Therapeutics", "Medicine", "Health", "Sciences", "Diagnostics"


