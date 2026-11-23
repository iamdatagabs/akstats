from akstats.ncbi import NCBIGEOClient


def test_ncbi_client_uses_optional_environment_credentials_without_hardcoding():
    client = NCBIGEOClient(email="research@example.org", api_key="test-key")
    params = client._eutils_params({"db": "gds", "term": "GSE158634[ACCN]"})
    assert params["tool"] == "akstats"
    assert params["email"] == "research@example.org"
    assert params["api_key"] == "test-key"
    assert params["db"] == "gds"


def test_ncbi_client_does_not_require_an_api_key_for_standard_access(monkeypatch):
    monkeypatch.delenv("NCBI_API_KEY", raising=False)
    client = NCBIGEOClient(api_key="")
    params = client._eutils_params({"db": "gds"})
    assert "api_key" not in params
