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



def test_ncbi_client_fetches_biosample_attributes_as_xml():
    class Response:
        content = (
            b'<BioSampleSet><BioSample accession="SAMN00000001">'
            b'<Description>SKH-1 sample</Description><Attributes>'
            b'<Attribute attribute_name="animal_id">Mouse 01</Attribute>'
            b'<Attribute attribute_name="treatment">UVR</Attribute>'
            b'</Attributes></BioSample></BioSampleSet>'
        )

        def raise_for_status(self):
            return None

    class Session:
        def __init__(self):
            self.headers = {}

        def mount(self, *_args):
            return None

        def get(self, *_args, **_kwargs):
            return Response()

    client = NCBIGEOClient(session=Session())
    records = client.fetch_biosample_records(["SAMN00000001"])
    assert records["SAMN00000001"]["title"] == "SKH-1 sample"
    assert records["SAMN00000001"]["attributes"] == {
        "animal_id": ["Mouse 01"],
        "treatment": ["UVR"],
    }



def test_ncbi_client_fetches_sra_instrument_and_biosample_link():
    class Response:
        content = (
            b'<EXPERIMENT_PACKAGE_SET><EXPERIMENT_PACKAGE>'
            b'<EXPERIMENT accession="SRX0000001"><PLATFORM><ILLUMINA>'
            b'<INSTRUMENT_MODEL>Illumina HiSeq 2500</INSTRUMENT_MODEL>'
            b'</ILLUMINA></PLATFORM></EXPERIMENT>'
            b'<SAMPLE accession="SRS00000001"><IDENTIFIERS>'
            b'<EXTERNAL_ID namespace="BioSample">SAMN00000001</EXTERNAL_ID>'
            b'</IDENTIFIERS></SAMPLE></EXPERIMENT_PACKAGE>'
            b'</EXPERIMENT_PACKAGE_SET>'
        )

        def raise_for_status(self):
            return None

    class Session:
        def __init__(self):
            self.headers = {}

        def mount(self, *_args):
            return None

        def get(self, *_args, **_kwargs):
            return Response()

    client = NCBIGEOClient(session=Session())
    records = client.fetch_sra_experiments(["SRX0000001"])
    assert records["SRX0000001"]["instrument_model"] == "Illumina HiSeq 2500"
    assert records["SRX0000001"]["sra_sample_accession"] == "SRS00000001"
    assert records["SRX0000001"]["biosample_accession"] == "SAMN00000001"



def test_ncbi_client_paces_entrez_requests_to_documented_limits(monkeypatch):
    monkeypatch.delenv("NCBI_API_KEY", raising=False)
    standard = NCBIGEOClient(api_key="")
    keyed = NCBIGEOClient(api_key="test-key")
    assert standard._eutils_interval == 1 / 3
    assert keyed._eutils_interval == 1 / 10
