"""pytest configuration for the raredx web tests."""

def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "network: test hits live REST APIs (Ensembl/gnomAD/ClinVar); needs internet",
    )
