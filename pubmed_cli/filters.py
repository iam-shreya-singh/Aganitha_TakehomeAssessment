def is_non_academic_journal(affiliation: str) -> bool:

academic_keywords = ["university","college", "institute","hospital", "school", "center"]
company_keywords =["ltd", "inc", "corp", "pharma", "biotech", "biotechnology", "pharmaceutical", "laboratory"]


aff_lower = affiliation.lower()
if any(keyword in aff_lower for keyword in academic_keywords):
    return False
if any(keyword in aff_lower for keyword in company_keywords):
    return True
return False

