# pubmed_fetcher/fetcher.py
import re
from typing import List, Dict, Optional, Set
from Bio import Entrez
import xml.etree.ElementTree as ET
import csv
import logging

logger = logging.getLogger(__name__)

ACADEMIC_KEYWORDS: Set[str] = {
    "university", "college", "institute", "school", "hospital", "dept", 
    "department", "research", "academy", "laboratory", "lab", "center", 
    "institut", "universidad", "université", "universität", "universitat", 
    "universitaet", "universidade", "universidad", "universitet", "universiteit", 
    "universitetet", "universiteti", "universitetet", "universiteti", 
    "universitetet", "universiteti", "universitetet", "universiteti"
}

PHARMA_KEYWORDS: Set[str] = {
    "pharma", "pharmaceutical", "biotech", "biotechnology", "inc", "llc", 
    "ltd", "sa", "ag", "co", "corp", " Holding", "Holdings", " GmbH", 
    " AB", " S.A.", "Pharmaceuticals", "Biotechnology", "Bioscience", 
    "Therapeutics", "Medicine", "Health", "Sciences", "Diagnostics"
}

class PubMedFetcher:
    def __init__(self, email: str, debug: bool = False):
        self.email = email
        self.debug = debug
        Entrez.email = email
        if debug:
            logging.basicConfig(level=logging.DEBUG)
        else:
            logging.basicConfig(level=logging.INFO)

    def fetch_papers(self, query: str, max_results: int = 100) -> List[Dict]:
        """Fetch papers from PubMed based on query"""
        logger.info(f"Fetching papers for query: {query}")
        
        try:
            handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
            record = Entrez.read(handle)
            ids = record["IdList"]
            
            if not ids:
                logger.info("No papers found")
                return []
                
            logger.info(f"Found {len(ids)} papers")
            handle = Entrez.efetch(db="pubmed", id=ids, retmode="xml")
            xml_data = handle.read()
            
            root = ET.fromstring(xml_data)
            papers = []
            
            for article in root.findall(".//PubmedArticle"):
                try:
                    paper_dict = self.process_paper(article)
                    if paper_dict:
                        papers.append(paper_dict)
                except Exception as e:
                    logger.error(f"Error processing paper: {e}")
                    if self.debug:
                        logger.exception("Full traceback:")
                    
            return papers
            
        except Exception as e:
            logger.error(f"Error fetching papers: {e}")
            if self.debug:
                logger.exception("Full traceback:")
            return []

    def process_paper(self, article) -> Optional[Dict]:
        """Process a single PubMed article"""
        try:
            pmid_elem = article.find(".//PMID")
            pmid = pmid_elem.text if pmid_elem is not None else ""
            
            title_elem = article.find(".//ArticleTitle")
            title = title_elem.text if title_elem is not None else ""
            
            pub_date = self._extract_publication_date(article)
            
            authors = article.findall(".//Author")
            non_academic_authors = set()
            company_affiliations = set()
            
            for author in authors:
                last_name = author.findtext("LastName", "")
                fore_name = author.findtext("ForeName", "")
                full_name = f"{fore_name} {last_name}".strip()
                
                if not full_name:
                    continue
                    
                # Check affiliations for this author
                affiliations = []
                aff_elem = author.find("Affiliation")
                if aff_elem is not None:
                    affiliations.append(aff_elem.text)
                
                for aff in affiliations:
                    if self._is_academic_affiliation(aff):
                        continue
                    if self._is_pharma_affiliation(aff):
                        non_academic_authors.add(full_name)
                        company_affiliations.add(aff.strip())
            
            # Find corresponding author email
            corresponding_email = self._find_corresponding_email(article, authors)
            
            return {
                "PubmedID": pmid,
                "Title": title,
                "Publication Date": pub_date,
                "Non-academic Author(s)": ", ".join(sorted(non_academic_authors)),
                "Company Affiliation(s)": ", ".join(sorted(company_affiliations)),
                "Corresponding Author Email": corresponding_email or ""
            }
            
        except Exception as e:
            logger.error(f"Error processing paper: {e}")
            if self.debug:
                logger.exception("Full traceback:")
            return None

    def _extract_publication_date(self, article) -> str:
        """Extract publication date in YYYY-MM-DD format"""
        pub_date_elem = article.find(".//PubDate")
        if pub_date_elem is None:
            return ""
            
        year_elem = pub_date_elem.find("Year")
        month_elem = pub_date_elem.find("Month")
        day_elem = pub_date_elem.find("Day")
        
        year = year_elem.text if year_elem is not None else ""
        month = month_elem.text if month_elem is not None else "01"
        day = day_elem.text if day_elem is not None else "01"
        
        # Handle month abbreviations
        month_map = {
            "Jan": "01", "Feb": "02", "Mar": "03", "Apr": "04",
            "May": "05", "Jun": "06", "Jul": "07", "Aug": "08",
            "Sep": "09", "Oct": "10", "Nov": "11", "Dec": "12"
        }
        
        if month in month_map:
            month = month_map[month]
        elif month.isdigit() and len(month) == 1:
            month = f"0{month}"
            
        if day.isdigit() and len(day) == 1:
            day = f"0{day}"
            
        return f"{year}-{month}-{day}"

    def _is_academic_affiliation(self, affiliation: str) -> bool:
        """Check if affiliation is academic"""
        if not affiliation:
            return False
            
        aff_lower = affiliation.lower()
        return any(keyword in aff_lower for keyword in ACADEMIC_KEYWORDS)

    def _is_pharma_affiliation(self, affiliation: str) -> bool:
        """Check if affiliation is pharma/biotech"""
        if not affiliation:
            return False
            
        aff_lower = affiliation.lower()
        return any(keyword in aff_lower for keyword in PHARMA_KEYWORDS)

    def _find_corresponding_email(self, article, authors) -> str:
        """Find corresponding author email"""
        # First try to find author marked as corresponding
        for author in authors:
            aff_elem = author.find("Affiliation")
            if aff_elem is not None and "corresponding" in aff_elem.text.lower():
                email_elem = author.find("Email")
                if email_elem is not None:
                    return email_elem.text
        
        # If not found, return first author email
        for author in authors:
            email_elem = author.find("Email")
            if email_elem is not None:
                return email_elem.text
                
        return ""

    def write_csv(self, papers: List[Dict], output_file: Optional[str] = None):
        """Write papers to CSV file or print to console"""
        if not papers:
            logger.warning("No papers to write")
            return
            
        fieldnames = [
            "PubmedID", "Title", "Publication Date", 
            "Non-academic Author(s)", "Company Affiliation(s)", 
            "Corresponding Author Email"
        ]
        
        if output_file:
            logger.info(f"Writing results to {output_file}")
            with open(output_file, 'w', newline='', encoding='utf-8') as csvfile:
                writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                writer.writeheader()
                for paper in papers:
                    if paper:  # Skip None values
                        writer.writerow(paper)
        else:
            logger.info("Writing results to console")
            writer = csv.DictWriter(open('stdout', 'w', encoding='utf-8'), fieldnames=fieldnames)
            writer.writeheader()
            for paper in papers:
                if paper:
                    writer.writerow(paper)