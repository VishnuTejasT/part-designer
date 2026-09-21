"""Vercel Python serverless function: JSON API wrapping the plasmid_design
package for the web frontend (see index.html / vercel.json rewrite)."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from flask import Flask, abort, jsonify, request, send_from_directory  # noqa: E402

from plasmid_design.codon_usage import (  # noqa: E402
    CodonUsageError,
    DEFAULT_ORGANISM,
    ORGANISM_TABLES,
    load_codon_table,
)
from plasmid_design.design import design_cds  # noqa: E402
from plasmid_design.part_finder import PartFinderError, check_cds, find_parts  # noqa: E402
from plasmid_design.optimize import (  # noqa: E402
    DEFAULT_GC_BOUNDS,
    MAX_PROTEIN_LENGTH,
    Limits,
    OptimizationError,
    OptimizationRequest,
    StructuralRegion,
    optimize_cds,
)
from plasmid_design.optimize_report import build_optimization_report  # noqa: E402
from plasmid_design.report import build_report  # noqa: E402
from plasmid_design.reverse_translate import ProteinSequenceError  # noqa: E402

app = Flask(__name__, static_folder=None)  # our own /static route serves the repo-root static/ dir

PROJECT_ROOT = Path(__file__).resolve().parents[1]


@app.route("/", methods=["GET"])
def home():
    html = (PROJECT_ROOT / "index.html").read_text()
    return app.response_class(html, mimetype="text/html")


@app.route("/glossary", methods=["GET"])
def glossary():
    return app.response_class((PROJECT_ROOT / "glossary.html").read_text(), mimetype="text/html")


@app.route("/static/<path:name>", methods=["GET"])
def static_files(name):
    static_dir = PROJECT_ROOT / "static"
    if not (static_dir / name).is_file():
        abort(404)
    return send_from_directory(static_dir, name, max_age=300)


TEMPERATURE_RANGE = (4.0, 45.0)


@app.route("/api/limits", methods=["GET"])
def limits():
    """Values the client needs for validation and for its 'defaults' comparison
    (the client only sends settings that differ from these)."""
    d = Limits()
    return jsonify({
        "max_protein_length": MAX_PROTEIN_LENGTH,
        "temperature_range": list(TEMPERATURE_RANGE),
        "defaults": {
            "temperature": 37.0,
            "cai_floor": d.cai_floor,
            "rare_codon_cutoff": d.rare_codon_w,
            "gc_min": DEFAULT_GC_BOUNDS[0],
            "gc_max": DEFAULT_GC_BOUNDS[1],
            "longest_allowed_stem": d.max_stem_bp - 1,
            "worst_window_energy": d.max_window_mfe,
            "start_region_energy": d.init_dg_floor,
        },
    })


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

    def number(name, default, label):
        raw = payload.get(name)
        if raw in (None, ""):
            return default, None
        try:
            return float(raw), None
        except (TypeError, ValueError):
            return default, f"{label} must be a number"

    temperature_c, err = number("temperature", 37.0, "temperature")
    if err:
        return jsonify({"error": err}), 400
    if not TEMPERATURE_RANGE[0] <= temperature_c <= TEMPERATURE_RANGE[1]:
        return jsonify({"error": "Choose a temperature between 4 and 45 \u00b0C."}), 400

    d = Limits()
    values = {}
    for key, default in (
        ("cai_floor", d.cai_floor), ("rare_codon_cutoff", d.rare_codon_w),
        ("longest_allowed_stem", d.max_stem_bp - 1), ("worst_window_energy", d.max_window_mfe),
        ("start_region_energy", d.init_dg_floor), ("gc_min", DEFAULT_GC_BOUNDS[0]), ("gc_max", DEFAULT_GC_BOUNDS[1]),
    ):
        values[key], err = number(key, default, key)
        if err:
            return jsonify({"error": err}), 400
    if not 0 < values["cai_floor"] <= 1:
        return jsonify({"error": "cai_floor must be between 0 and 1"}), 400
    if not 0 <= values["rare_codon_cutoff"] <= 1:
        return jsonify({"error": "rare_codon_cutoff must be between 0 and 1"}), 400
    if not 0 <= values["gc_min"] < values["gc_max"] <= 1:
        return jsonify({"error": "gc_min must be lower than gc_max, both between 0 and 1"}), 400
    if values["longest_allowed_stem"] < 1:
        return jsonify({"error": "longest_allowed_stem must be at least 1"}), 400
    req_limits = Limits(
        cai_floor=values["cai_floor"], rare_codon_w=values["rare_codon_cutoff"],
        max_stem_bp=int(values["longest_allowed_stem"]) + 1,
        max_window_mfe=values["worst_window_energy"], init_dg_floor=values["start_region_energy"],
    )

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
            limits=req_limits,
            gc_bounds=(values["gc_min"], values["gc_max"]),
            vector_provides_start=bool(payload.get("vector_provides_start", False)),
        )
        results = optimize_cds(req)
    except (ProteinSequenceError, CodonUsageError, OptimizationError, ValueError) as exc:
        return jsonify({"error": str(exc)}), 400

    return jsonify(build_optimization_report(results, req_limits))


@app.route("/api/part-finder", methods=["POST"])
def part_finder():
    """Promoter / RBS / terminator suggestions from the verified Registry extract.
    Body: {host, standard?, level?, top?, kinds?, avoid_type_iis?, cds?}"""
    payload = request.get_json(silent=True) or {}
    host = payload.get("host")
    if not host:
        return jsonify({"error": "host is required"}), 400
    try:
        top = int(payload.get("top", 5))
    except (TypeError, ValueError):
        return jsonify({"error": "top must be a whole number"}), 400
    try:
        result = find_parts(
            host, standard=payload.get("standard") or "RFC10", level=payload.get("level") or None, top=top,
            kinds=tuple(payload.get("kinds") or ("promoter", "rbs", "terminator")),
            avoid_type_iis=bool(payload.get("avoid_type_iis", False)),
        )
    except PartFinderError as exc:
        return jsonify({"error": str(exc)}), 400
    cds = "".join(str(payload.get("cds") or "").split()).upper()
    if cds:
        if set(cds) - set("ACGT"):
            return jsonify({"error": "cds must contain only A, C, G, T"}), 400
        result["cds_rfc10"] = check_cds(cds)
    return jsonify(result)


@app.route("/api/health", methods=["GET"])
def health():
    # Touches the codon table loader so a broken data file surfaces here
    # instead of silently on first user request.
    load_codon_table(DEFAULT_ORGANISM)
    return jsonify({"status": "ok"})
