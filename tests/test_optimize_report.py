from plasmid_design.optimize import OptimizationRequest, optimize_cds
from plasmid_design.optimize_report import build_optimization_report, format_optimization_report_text

TEST_PROTEIN = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQ"


def _results():
    req = OptimizationRequest(
        protein=TEST_PROTEIN, hosts=("e_coli_k12", "synechocystis_pcc6803"), mode="PRODUCTION", seed=1
    )
    return optimize_cds(req)


def test_report_dict_structure():
    report = build_optimization_report(_results())
    assert len(report["sequences"]) == 2
    assert len(report["summary_table"]) == 2
    assert "e_coli_k12" in report["codon_table_sources"]
    assert report["trna_table_sources"]["synechocystis_pcc6803"] is None
    assert report["trna_table_sources"]["e_coli_k12"] is not None
    assert "in vitro" in report["caveat"]


def test_text_report_contains_caveat_and_sources():
    report = build_optimization_report(_results())
    text = format_optimization_report_text(report)
    assert "SUMMARY" in text
    assert report["caveat"] in text
    assert "Genomic tRNA Database" in text
    assert "Kazusa" in text
