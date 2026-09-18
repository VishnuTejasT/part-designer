"""Vercel Python serverless function: JSON API wrapping the plasmid_design
package for the web frontend (see index.html / vercel.json rewrite)."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from flask import Flask, jsonify, request  # noqa: E402

from plasmid_design.codon_usage import (  # noqa: E402
    CodonUsageError,
    DEFAULT_ORGANISM,
    ORGANISM_TABLES,
    load_codon_table,
)
from plasmid_design.design import design_cds  # noqa: E402
from plasmid_design.optimize import (  # noqa: E402
    OptimizationError,
    OptimizationRequest,
    StructuralRegion,
    optimize_cds,
)
from plasmid_design.optimize_report import build_optimization_report  # noqa: E402
from plasmid_design.report import build_report  # noqa: E402
from plasmid_design.reverse_translate import ProteinSequenceError  # noqa: E402

app = Flask(__name__)

PROJECT_ROOT = Path(__file__).resolve().parents[1]


@app.route("/", methods=["GET"])
def home():
    html = (PROJECT_ROOT / "index.html").read_text()
    return app.response_class(html, mimetype="text/html")


@app.route("/api/organisms", methods=["GET"])
def organisms():
    return jsonify({"organisms": sorted(ORGANISM_TABLES), "default": DEFAULT_ORGANISM})


@app.route("/api/convert", methods=["POST"])
def convert():
    payload = request.get_json(silent=True) or {}
    protein = payload.get("protein", "")
    organism = payload.get("organism") or DEFAULT_ORGANISM

    raw_seed = payload.get("seed")
    seed = None
    if raw_seed not in (None, ""):
        try:
            seed = int(raw_seed)
        except (TypeError, ValueError):
            return jsonify({"error": "seed must be an integer"}), 400

    if organism not in ORGANISM_TABLES:
        available = ", ".join(sorted(ORGANISM_TABLES))
        return jsonify({"error": f"unknown organism '{organism}'. Available: {available}"}), 400

    auto_fix = payload.get("auto_fix", True)

    try:
        result = design_cds(protein, organism=organism, seed=seed, auto_fix=auto_fix)
    except (ProteinSequenceError, CodonUsageError) as exc:
        return jsonify({"error": str(exc)}), 400

    report = build_report(result.translation, result.validation, result.repair)
    return jsonify(report)


@app.route("/api/optimize", methods=["POST"])
def optimize():
    payload = request.get_json(silent=True) or {}
    protein = payload.get("protein", "")
    hosts = payload.get("hosts") or []
    if isinstance(hosts, str):
        hosts = [h.strip() for h in hosts.split(",") if h.strip()]
    hosts = tuple(hosts)
    if not hosts:
        return jsonify({"error": "at least one host is required"}), 400

    mode = (payload.get("mode") or "BALANCED").upper()

    raw_seed = payload.get("seed")
    seed = None
    if raw_seed not in (None, ""):
        try:
            seed = int(raw_seed)
        except (TypeError, ValueError):
            return jsonify({"error": "seed must be an integer"}), 400

    forbidden_raw = payload.get("forbidden_enzymes") or []
    try:
        forbidden = tuple((item["name"], item["pattern"].upper()) for item in forbidden_raw)
    except (KeyError, TypeError, AttributeError):
        return jsonify({"error": "forbidden_enzymes must be a list of {name, pattern}"}), 400

    regions_raw = payload.get("structural_regions") or []
    try:
        regions = tuple(
            StructuralRegion(start=r["start"], end=r["end"], kind=r["kind"]) for r in regions_raw
        )
    except (KeyError, TypeError):
        return jsonify({"error": "structural_regions must be a list of {start, end, kind}"}), 400

    raw_temperature = payload.get("temperature")
    temperature_c = 37.0
    if raw_temperature not in (None, ""):
        try:
            temperature_c = float(raw_temperature)
        except (TypeError, ValueError):
            return jsonify({"error": "temperature must be a number"}), 400

    try:
        req = OptimizationRequest(
            protein=protein,
            hosts=hosts,
            mode=mode,
            native_organism=payload.get("native_organism") or None,
            five_prime_utr=payload.get("five_prime_utr", "") or "",
            n_tag=payload.get("n_tag", "") or "",
            c_tag=payload.get("c_tag", "") or "",
            forbidden_enzymes=forbidden,
            structural_regions=regions,
            hedge=bool(payload.get("hedge", False)),
            seed=seed,
            temperature_c=temperature_c,
        )
        results = optimize_cds(req)
    except (ProteinSequenceError, CodonUsageError, OptimizationError, ValueError) as exc:
        return jsonify({"error": str(exc)}), 400

    return jsonify(build_optimization_report(results))


@app.route("/api/health", methods=["GET"])
def health():
    # Touches the codon table loader so a broken data file surfaces here
    # instead of silently on first user request.
    load_codon_table(DEFAULT_ORGANISM)
    return jsonify({"status": "ok"})
