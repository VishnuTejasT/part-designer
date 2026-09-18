import json

from plasmid_design.cli import main


def test_cli_clean_sequence_exits_zero(capsys):
    exit_code = main(["--protein", "MKT", "--seed", "1"])
    assert exit_code == 0
    out = capsys.readouterr().out
    assert "RFC10 COMPATIBILITY VALIDATION" in out


def test_cli_exits_nonzero_on_invalid_protein(capsys):
    exit_code = main(["--protein", "MKT123", "--seed", "1"])
    assert exit_code == 2
    err = capsys.readouterr().err
    assert "error:" in err


def test_cli_json_out_written(tmp_path, capsys):
    out_path = tmp_path / "report.json"
    exit_code = main(
        ["--protein", "MKT", "--seed", "1", "--json-out", str(out_path), "--quiet"]
    )
    assert exit_code == 0
    captured = capsys.readouterr().out
    assert captured == ""  # --quiet suppresses stdout report

    with out_path.open() as fh:
        report = json.load(fh)
    assert report["protein"]["sequence"] == "MKT"
