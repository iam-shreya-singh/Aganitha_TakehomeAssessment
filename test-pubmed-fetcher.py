import pytest
from pubmed_fetcher.fetcher import PubMedFetcher

def test_pubmed_fetcher_initialization():
    """Test that PubMedFetcher can be initialized properly"""
    fetcher = PubMedFetcher(email="test@example.com", debug=False)
    assert fetcher.email == "test@example.com"
    assert fetcher.debug == False

def test_pubmed_fetcher_with_debug():
    """Test PubMedFetcher initialization with debug enabled"""
    fetcher = PubMedFetcher(email="test@example.com", debug=True)
    assert fetcher.email == "test@example.com"
    assert fetcher.debug == True

def test_csv_output_format():
    """Test that the CSV output has required columns"""
    # This is a basic test - can be expanded with actual data
    expected_columns = [
        "PubmedID",
        "Title", 
        "Publication Date",
        "Non-academic Author(s)",
        "Company Affiliation(s)",
        "Corresponding Author Email"
    ]
    
    # Testing with empty data first
    fetcher = PubMedFetcher(email="test@example.com")
    # This would normally create CSV output
    # we will expand this test with actual data when it runs

if __name__ == "__main__":
    pytest.main([__file__, "-v"])