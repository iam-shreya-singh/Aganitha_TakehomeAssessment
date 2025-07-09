from Bio import Entrez
from xml.etree import ElementTree as ET

Entrez.email = "shreyasinghin24@gmail.com"

def fetch_pubmed(query, max_results=10):
    """
    Fetch PubMed articles based on user input query...
    """
    handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
    record = Entrez.read(handle)
    ids = record.get("IdList", [])
    if debug:
        print(f"[DEBUG] Found PubMed ID's: {ids}")
    
    if not ids:
        print("No articles/ papers found for the given query.")
        return []
    
    handle = Entrez.efetch(db="pubmed", id=",".join(ids), retmode="xml")
    xml_data = handle.read()
    root = ET.fromstring(xml_data)

    papers = []
    for article in root.findall(".//PubmedArticle"):
       pmid = article.findtext(".//PMID")
       title = aricle.findtext(".//ArticleTitle")
       pub_date = article.findtext(".//PubDate/Year") or "Unknown"
       authors = [
           f" {a.findtext('ForeName', '')} {a.findtext('LastName', '')}".strip()
           for a in article.findall(".//Author") if a.findtext("LastName")

  ]
       affiliation = [
           aff.text
           for aff in article.findall(".//AffiliationInfo/Affiliation")
           if aff.text
       ]

       papers.append({
           "PubMed ID": pmid,
           "Title": title,
           "PublicationDate": pub_date,
           "Authors": authors,
           "Affiliation": affiliation
       })
       return papers
    