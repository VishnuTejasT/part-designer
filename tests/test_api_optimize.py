import json

import pytest

flask = pytest.importorskip("flask")

from api.index import app  # noqa: E402

UBQ = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"


@pytest.fixture()
def client():
    return app.test_client()


def post(client, **body):
    return client.post("/api/optimize", json={"protein": UBQ, "hosts": ["e_coli_k12"], "mode": "PRODUCTION",
                                              "seed": 1, **body})


def test_limits_endpoint(client):
    d = client.get("/api/limits").get_json()
    assert d["max_protein_length"] == 1500
    assert d["temperature_range"] == [4.0, 45.0]
    assert d["defaults"]["cai_floor"] == 0.9 and d["defaults"]["longest_allowed_stem"] == 7


def test_default_request_matches_explicit_defaults(client):
    a = post(client).get_json()["sequences"][0]["dna"]
    b = post(client, cai_floor=0.9, rare_codon_cutoff=0.3, gc_min=0.3, gc_max=0.7,
             longest_allowed_stem=7, worst_window_energy=-15, start_region_energy=-5).get_json()["sequences"][0]["dna"]
    assert a == b


def test_report_has_checks_and_limits(client):
    d = post(client, cai_floor=0.85).get_json()
    assert d["limits_used"]["cai_floor"] == 0.85
    assert d["sequences"][0]["checks_total"] == 10


@pytest.mark.parametrize("temp", [3.9, 45.1, 60])
def test_temperature_out_of_range(client, temp):
    r = post(client, temperature=temp)
    assert r.status_code == 400 and "between 4 and 45" in r.get_json()["error"]


def test_too_long_protein_message(client):
    r = post(client, protein="M" + "A" * 1500)
    assert r.status_code == 400
    assert r.get_json()["error"].endswith("The limit is 1,500. Try splitting it into domains.")


def test_vector_provides_start(client):
    plain = post(client).get_json()["sequences"][0]["dna"]
    trimmed = post(client, vector_provides_start=True).get_json()["sequences"][0]["dna"]
    assert trimmed == plain[3:]


def test_bad_limit_values(client):
    assert post(client, cai_floor="abc").status_code == 400
    assert post(client, gc_min=0.8, gc_max=0.2).status_code == 400


def test_static_missing_is_404(client):
    assert client.get("/static/nope.js").status_code == 404


@pytest.mark.parametrize("path,mime", [("/static/app.js", "javascript"), ("/static/logic.js", "javascript"),
                                       ("/static/strings.js", "javascript"), ("/static/app.css", "css")])
def test_static_assets_are_served(client, path, mime):
    r = client.get(path)
    assert r.status_code == 200 and mime in r.mimetype and len(r.data) > 100


def test_pages_are_served(client):
    assert b"/static/app.js" in client.get("/").data
    assert b"PD_STRINGS" in client.get("/glossary").data


def test_part_finder_endpoint(client):
    r = client.post("/api/part-finder", json={"host": "e_coli_bl21_de3", "level": "high", "top": 2, "cds": "ATGGAATTCTAA"})
    d = r.get_json()
    assert r.status_code == 200 and set(d["results"]) == {"promoter", "rbs", "terminator"}
    assert all(len(b["parts"]) <= 2 for b in d["results"].values())
    assert d["cds_rfc10"]["illegal_sites"] == ["EcoRI"]
    assert d["dataset"]["source"]["sha256"].startswith("c64bbc9a")


@pytest.mark.parametrize("body", [{}, {"host": "c_reinhardtii"}, {"host": "e_coli_k12", "top": "x"},
                                  {"host": "e_coli_k12", "level": "extreme"}, {"host": "e_coli_k12", "cds": "ATGXX"}])
def test_part_finder_rejects_bad_input(client, body):
    assert client.post("/api/part-finder", json=body).status_code == 400
