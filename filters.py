def is_non_academic(affiliation: str) -> bool:
    """
    Determine if an affiliation is non-academic.

    Args:
        affiliation (str): The affiliation string from PubMed.

    Returns:
        bool: True if the affiliation is non-academic (company), False if academic.
    """
    academic_keywords = [
        "university", "college", "institute", "hospital", "school", "center"
    ]
    company_keywords = [
        "ltd", "inc", "corp", "pharma", "biotech", "biotechnology", "pharmaceutical", "laboratory"
    ]

    aff_lower = affiliation.lower()

    # Check if affiliation contains academic keywords
    if any(word in aff_lower for word in academic_keywords):
        return False

    # Check if affiliation contains company keywords
    if any(word in aff_lower for word in company_keywords):
        return True

    # Default: treat as academic unless clear company keywords
    return False
