from Bio import Entrez
from xml.etree import ElementTree as ET

Entrez.email = "shreyasinghin24@gmail.com"

def fetch_pubmed(query, max_results=10, debug=False):
    """
    Fetch PubMed articles based on user input query.

    Args:
        query (str): Search term.
        max_results (int): Max number of papers to fetch.
        debug (bool): Enable debug output.

    Returns:
        list: A list of dictionaries containing paper details.
    """
    if debug:
        print(f"[DEBUG] Searching PubMed for query: {query}")

    # Search for paper IDs
    handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
    record = Entrez.read(handle)
    ids = record.get("IdList", [])
    handle.close()

    if debug:
        print(f"[DEBUG] Found PubMed IDs: {ids}")

    if not ids:
        print("No articles/papers found for the given query.")
        return []

    # Fetch paper details
    handle = Entrez.efetch(db="pubmed", id=",".join(ids), retmode="xml")
    xml_data = handle.read()
    handle.close()
    root = ET.fromstring(xml_data)

    papers = []
    for article in root.findall(".//PubmedArticle"):
        pmid = article.findtext(".//PMID")
        title = article.findtext(".//ArticleTitle")
        pub_date = article.findtext(".//PubDate/Year") or "Unknown"
        authors = [
            f"{a.findtext('ForeName', '')} {a.findtext('LastName', '')}".strip()
            for a in article.findall(".//Author") if a.findtext("LastName")
        ]
        affiliations = [
            aff.text
            for aff in article.findall(".//AffiliationInfo/Affiliation")
            if aff.text
        ]

        papers.append({
            "PubmedID": pmid,
            "Title": title,
            "PublicationDate": pub_date,
            "Authors": authors,
            "Affiliations": affiliations
        })

    if debug:
        print(f"[DEBUG] Retrieved {len(papers)} papers.")
    return papers
