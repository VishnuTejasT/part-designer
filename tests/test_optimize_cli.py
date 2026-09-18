import json

from plasmid_design.optimize_cli import main

TEST_PROTEIN = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQ"


def test_cli_runs_and_passes(capsys):
    exit_code = main(["--protein", TEST_PROTEIN, "--hosts", "e_coli_k12", "--mode", "PRODUCTION", "--seed", "1"])
    out = capsys.readouterr().out
    assert exit_code == 0
    assert "HOST: e_coli_k12" in out
    assert "SUMMARY" in out


def test_cli_json_out(tmp_path):
    out_path = tmp_path / "report.json"
    exit_code = main(
        [
            "--protein", TEST_PROTEIN,
            "--hosts", "e_coli_k12,synechocystis_pcc6803",
            "--mode", "PRODUCTION",
            "--seed", "2",
            "--quiet",
            "--json-out", str(out_path),
        ]
    )
    assert exit_code == 0
    data = json.loads(out_path.read_text())
    assert len(data["summary_table"]) == 2


def test_cli_invalid_protein_returns_error_code(capsys):
    exit_code = main(["--protein", "not a protein 123", "--hosts", "e_coli_k12"])
    err = capsys.readouterr().err
    assert exit_code == 2
    assert "error:" in err


def test_cli_forbidden_enzymes_flag():
    exit_code = main(
        [
            "--protein", TEST_PROTEIN,
            "--hosts", "e_coli_k12",
            "--mode", "PRODUCTION",
            "--forbidden-enzymes", "EcoRV:GATATC",
            "--seed", "3",
            "--quiet",
        ]
    )
    assert exit_code == 0
