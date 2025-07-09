def is_non_academic(affiliation: str) -> bool:

academic_keywords = ["university","college", "institute","hospital", "school", "center"]
company_keywords =["ltd", "inc", "corp", "pharma", "biotech", "biotechnology", "pharmaceutical", "laboratory"]


aff_lower = affiliation.lower()
if any(word in aff_lower for word in academic_keywords):
    return False
if any(word in aff_lower for word in company_keywords):
    return True
return False

