/* Pure client logic: sequence cleaning and validation, request building,
   result merging. No DOM access, so node can test it. All display text comes
   from strings.js. */
(function (root) {
  "use strict";
  var S = typeof require === "function" && typeof module !== "undefined" ? require("./strings.js") : root.PD_STRINGS;

  var AMINO = "ACDEFGHIKLMNPQRSTVWY";
  var AMBIGUOUS = { B: "kindAmbiguous", Z: "kindAmbiguous", J: "kindAmbiguous", X: "kindX", U: "kindRare", O: "kindRare" };
  var IUPAC = "ACGTRYSWKMBDHVN";
  var CODONS = (function () {
    var bases = "TCAG", aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", t = {}, i = 0;
    for (var a = 0; a < 4; a++) for (var b = 0; b < 4; b++) for (var c = 0; c < 4; c++) t[bases[a] + bases[b] + bases[c]] = aas[i++];
    return t;
  })();

  function plural(n) { return n.toLocaleString("en-US"); }

  /* ---- FASTA / cleaning --------------------------------------------------- */

  function splitRecords(raw) {
    var lines = String(raw).replace(/\r/g, "").split("\n");
    var records = [], cur = null;
    lines.forEach(function (line) {
      if (line.trim().charAt(0) === ">") {
        cur = { name: line.trim().slice(1).trim(), lines: [] };
        records.push(cur);
      } else {
        if (!cur) { cur = { name: "", lines: [], noHeader: true }; records.push(cur); }
        cur.lines.push(line);
      }
    });
    return records.filter(function (r) { return r.lines.join("").trim() !== "" || r.name; });
  }

  /* Clean one pasted block. Returns the sequence plus everything the UI needs
     to explain what happened. `chosenRecord` picks a record from a multi-FASTA. */
  function analyzeProtein(raw, opts) {
    opts = opts || {};
    var limit = opts.maxLength || 1500;
    var out = { sequence: "", name: "", notes: [], errors: [], warnings: [], infos: [], invalid: [],
                records: 0, dna: false, blocked: false, needsRecordChoice: false };
    raw = String(raw == null ? "" : raw);
    if (raw.trim() === "") { out.blocked = true; out.empty = true; return out; }

    var records = splitRecords(raw);
    out.records = records.length;
    if (records.length > 1 && opts.chosenRecord === undefined) {
      out.needsRecordChoice = true; out.blocked = true;
      out.errors.push({ id: "multi", text: S.fmt(S.errors.multi, { n: records.length }) });
      out.recordNames = records.map(function (r, i) { return r.name || ("#" + (i + 1)); });
      return out;
    }
    var rec = records[Math.min(opts.chosenRecord || 0, records.length - 1)] || { name: "", lines: [""] };
    if (rec.name) { out.name = rec.name; out.notes.push({ id: "fasta", text: S.errors.fastaHeader }); }

    var body = rec.lines.join(rec.noHeader && rec.lines.length > 1 ? "\n" : "");
    var removedNoise = /[0-9 \t]/.test(body) || (rec.noHeader && /\n/.test(body.trim()));
    body = body.replace(/[0-9\s]/g, "").toUpperCase();
    if (removedNoise) out.notes.push({ id: "cleaned", text: S.errors.cleaned });

    if (/\*$/.test(body)) { body = body.replace(/\*+$/, ""); out.notes.push({ id: "star", text: S.errors.trailingStar }); }

    // DNA pasted instead of protein: >= 95% A, C, G, T or U.
    var letters = body.replace(/[^A-Z]/g, "");
    if (letters.length >= 9) {
      var dnaCount = (letters.match(/[ACGTU]/g) || []).length;
      if (dnaCount / letters.length >= 0.95) {
        out.dna = true; out.blocked = true; out.dnaSequence = letters;
        out.errors.push({ id: "dna", text: S.errors.dna });
        out.sequence = body;
        return out;
      }
    }

    if (body.indexOf("*") !== -1) {
      out.blocked = true;
      out.errors.push({ id: "middleStar", text: S.errors.middleStar });
      var at = body.indexOf("*");
      out.invalid.push({ index: at, char: "*" });
    }

    var perLetter = {}, generic = [];
    for (var i = 0; i < body.length; i++) {
      var ch = body[i];
      if (ch === "*" || AMINO.indexOf(ch) !== -1) continue;
      if (AMBIGUOUS[ch]) { perLetter[ch] = true; out.invalid.push({ index: i, char: ch }); }
      else { generic.push(i); out.invalid.push({ index: i, char: ch }); }
    }
    Object.keys(perLetter).forEach(function (l) {
      out.blocked = true;
      out.errors.push({ id: "badLetter", letter: l, text: S.fmt(S.errors.badLetterX, { letter: l, kind: S.errors[AMBIGUOUS[l]] }) });
    });
    if (generic.length) {
      out.blocked = true;
      out.errors.push({ id: "invalidLetters", removable: true,
        text: S.fmt(generic.length === 1 ? S.errors.invalidLettersOne : S.errors.invalidLetters, { n: generic.length }) });
    }

    out.sequence = body;
    var n = body.replace(/\*/g, "").length;
    if (n === 0) { out.blocked = true; out.empty = true; return out; }
    if (n < 3) { out.blocked = true; out.errors.push({ id: "short", text: S.fmt(S.errors.tooShortBlock, { n: n }) }); }
    else if (n < 20) out.warnings.push({ id: "short", text: S.fmt(S.errors.tooShortWarn, { n: n }) });
    if (n > limit) { out.blocked = true; out.errors.push({ id: "long", text: S.fmt(S.errors.tooLong, { n: plural(n), limit: plural(limit) }) }); }
    if (body[0] !== "M") out.infos.push({ id: "noMet", text: S.errors.noMet });
    return out;
  }

  /* Remove every character that is not an amino acid (the "Remove them" button). */
  function removeInvalid(sequence) {
    return sequence.split("").filter(function (c) { return AMINO.indexOf(c) !== -1; }).join("");
  }

  /* Reading frame 1, DNA or RNA, stopping at the first stop codon. */
  function translateFrame1(dna) {
    var s = String(dna).toUpperCase().replace(/U/g, "T").replace(/[^ACGT]/g, ""), aa = "";
    for (var i = 0; i + 3 <= s.length; i += 3) {
      var a = CODONS[s.substr(i, 3)];
      if (a === "*") break;
      aa += a;
    }
    return aa;
  }

  /* ---- Field validation ---------------------------------------------------- */

  function validateTag(aa) {
    aa = String(aa || "").toUpperCase().replace(/\s/g, "");
    if (!aa) return { ok: true, sequence: "" };
    for (var i = 0; i < aa.length; i++) {
      if (AMINO.indexOf(aa[i]) === -1) return { ok: false, sequence: aa, error: S.fmt(S.errors.tagBad, { tag: aa, letter: aa[i] }) };
    }
    return { ok: true, sequence: aa, warning: aa.length > 30 ? S.fmt(S.errors.tagLong, { n: aa.length }) : null };
  }

  function validateSite(site) {
    site = String(site || "").toUpperCase().replace(/\s/g, "");
    for (var i = 0; i < site.length; i++) {
      if (IUPAC.indexOf(site[i]) === -1) return { ok: false, error: S.fmt(S.errors.siteBad, { site: site, letter: site[i] }) };
    }
    if (site.length < 4 || site.length > 10) return { ok: false, error: S.errors.siteLength };
    return { ok: true, site: site };
  }

  function validateStartDna(text) {
    var t = String(text || "").toUpperCase().replace(/\s/g, "");
    if (!t) return { ok: true, dna: "" };
    if (/[^ACGTU]/.test(t)) return { ok: false, error: S.errors.startDna };
    var hadU = t.indexOf("U") !== -1;
    return { ok: true, dna: t.replace(/U/g, "T"), note: hadU ? S.errors.startDnaU : null };
  }

  function validateTemperature(v) {
    var n = Number(v);
    if (v === "" || v === null || isNaN(n) || n < 4 || n > 45) return { ok: false, error: S.errors.temperature };
    return { ok: true, value: n };
  }

  function validateSeed(v) {
    if (v === "" || v === null || v === undefined) return { ok: true, value: null };
    if (!/^-?\d+$/.test(String(v).trim())) return { ok: false, error: S.errors.seed };
    return { ok: true, value: parseInt(v, 10) };
  }

  /* Regions are 1-based, matching the engine. Returns per-row errors. */
  function validateRegions(regions, proteinLength) {
    var errors = regions.map(function () { return null; });
    regions.forEach(function (r, i) {
      if (!(r.start >= 1)) errors[i] = S.errors.regionStart;
      else if (r.end < r.start) errors[i] = S.errors.regionOrder;
      else if (r.end > proteinLength) errors[i] = S.fmt(S.errors.regionRange, { end: r.end, n: proteinLength });
    });
    for (var a = 0; a < regions.length; a++) for (var b = a + 1; b < regions.length; b++) {
      if (regions[a].start <= regions[b].end && regions[b].start <= regions[a].end) {
        if (!errors[a]) errors[a] = S.errors.regionOverlap;
        if (!errors[b]) errors[b] = S.errors.regionOverlap;
      }
    }
    return errors;
  }

  /* Likely flexible linkers: runs of 5+ glycine/serine (e.g. GGGGS repeats). */
  function detectRegions(protein) {
    var out = [], re = /[GS]{5,}/g, m;
    while ((m = re.exec(protein))) out.push({ start: m.index + 1, end: m.index + m[0].length, type: "linker", auto: true });
    return out;
  }

  /* ---- Requests ------------------------------------------------------------- */

  var MODE_FOR_GOAL = { PRODUCTION: "PRODUCTION", FOLDING: "FOLDING", BALANCED: "BALANCED", ALL: "ALL" };
  var DEFAULT_TEMPERATURE = 37;

  function hostTemperature(hostId) { return (S.hosts[hostId] || {}).temp || DEFAULT_TEMPERATURE; }

  /* Build the API requests for the current form state. The client sends only
     values the user set or that differ from engine defaults. Hosts with
     different suggested temperatures (and no user-set temperature) go in
     separate requests, because the API takes one temperature per request.
     `state`: see app.js. `defaults`: from /api/limits. */
  function buildRequests(state, defaults) {
    var protein = state.protein;
    if (protein[0] !== "M") protein = "M" + protein; // we add the start signal (or the vector's ATG is virtual)
    var base = { protein: protein, mode: MODE_FOR_GOAL[state.goal] };
    if (state.goal !== "PRODUCTION" && state.naturalOrganism) base.native_organism = state.naturalOrganism;
    if (state.startDna) base.five_prime_utr = state.startDna;
    if (state.nTag) base.n_tag = state.nTag;
    if (state.cTag) base.c_tag = state.cTag;
    if (state.vectorAtg) base.vector_provides_start = true;
    if (state.sites && state.sites.length) base.forbidden_enzymes = state.sites.map(function (s) { return { name: s.name, pattern: s.site }; });
    if (state.goal !== "PRODUCTION" && state.regions && state.regions.length) {
      base.structural_regions = state.regions.map(function (r) { return { start: r.start, end: r.end, kind: r.type }; });
    }
    if (state.hedge) base.hedge = true;
    if (state.seed !== null && state.seed !== undefined) base.seed = state.seed;

    var d = defaults || {};
    var lim = state.limits || {};
    [["cai_floor", "cai"], ["worst_window_energy", "mfe"], ["gc_min", "gcMin"], ["gc_max", "gcMax"], ["longest_allowed_stem", "stem"]]
      .forEach(function (p) {
        var v = lim[p[1]];
        if (v !== undefined && v !== null && v !== "" && Number(v) !== d[p[0]]) base[p[0]] = Number(v);
      });
    ["cai_floor", "worst_window_energy", "start_region_energy", "rare_codon_cutoff", "longest_allowed_stem"].forEach(function (k) {
      if (state.overrides && state.overrides[k] !== undefined) base[k] = state.overrides[k];
    });

    var groups = {};
    state.hosts.forEach(function (h) {
      var t = state.temperatureEdited ? Number(state.temperature) : hostTemperature(h);
      (groups[t] = groups[t] || []).push(h);
    });
    return Object.keys(groups).map(function (t) {
      var req = Object.assign({}, base, { hosts: groups[t] });
      if (Number(t) !== DEFAULT_TEMPERATURE) req.temperature = Number(t);
      return req;
    });
  }

  /* Combine several API reports into one (one per temperature group). */
  function mergeReports(reports) {
    if (reports.length === 1) return reports[0];
    var out = { sequences: [], summary_table: [], codon_table_sources: {}, trna_table_sources: {}, caveat: reports[0].caveat, limits_used: reports[0].limits_used };
    reports.forEach(function (r) {
      out.sequences = out.sequences.concat(r.sequences);
      out.summary_table = out.summary_table.concat(r.summary_table);
      Object.assign(out.codon_table_sources, r.codon_table_sources);
      Object.assign(out.trna_table_sources, r.trna_table_sources);
    });
    return out;
  }

  /* ---- Results helpers -------------------------------------------------------- */

  function statusCounts(seq) {
    var review = 0, fail = 0, pass = 0;
    seq.checks.forEach(function (c) {
      if (!c.applicable) return;
      if (c.status === "pass") pass++; else if (c.status === "review") review++; else fail++;
    });
    return { pass: pass, review: review, fail: fail, total: pass + review + fail };
  }

  function bannerText(seq) {
    var c = statusCounts(seq), open = c.review + c.fail;
    if (open === 0) return { ok: true, text: S.fmt(S.results.ready, { n: c.total }) };
    return { ok: false, text: open === 1 ? S.results.readyOne : S.fmt(S.results.readySome, { k: open }) };
  }

  function checkResultText(check) {
    var tpl = S.checks[check.id], v = check.value || {};
    var vars = {
      matched: v.matched, total: v.total, found: v.found, lowest: v.lowest !== undefined ? Number(v.lowest).toFixed(2) : "",
      cai: v.cai !== undefined && v.cai !== null ? Number(v.cai).toFixed(2) : "", goal: v.goal !== undefined ? Number(v.goal).toFixed(2) : "",
      // One decimal: rounding to a whole number made 29.6% read as a passing "30%" against a 30% limit.
      gcMin: v.gc_min !== undefined && v.gc_min !== null ? (v.gc_min * 100).toFixed(1) : "",
      gcMax: v.gc_max !== undefined && v.gc_max !== null ? (v.gc_max * 100).toFixed(1) : "",
      dG: v.dG !== undefined && v.dG !== null ? Number(v.dG).toFixed(1) : "",
      stem: v.longest_stem, strongest: v.strongest !== undefined && v.strongest !== null ? Number(v.strongest).toFixed(1) : "",
    };
    var key = check.status === "review" && tpl.review ? "review" : check.status === "pass" ? "pass" : "fail";
    return S.fmt(tpl[key], vars);
  }

  function techName(check) {
    var v = check.value || {};
    return S.fmt(S.checks[check.id].tech, { cutoff: v.cutoff !== undefined ? v.cutoff : 0.3 });
  }

  /* Turn a structured engine conflict into fix-panel content, or null when
     there is nothing the user can act on. */
  function fixPanel(seq, limitsUsed) {
    var c = (seq.conflicts || []).filter(function (x) { return x.constraint !== "sequence_rules"; })[0];
    if (!c) {
      // Repeats, GC balance and similar: the engine couldn't clear them, and only a new random try can.
      var rules = (seq.conflicts || []).filter(function (x) { return x.constraint === "sequence_rules"; })[0];
      if (!rules) return null;
      var plain = (rules.reasons || []).map(function (x) { return S.results.synthReasons[x] || S.results.synthOther; })
        .filter(function (x, i, a) { return a.indexOf(x) === i; });
      return { constraint: "sequence_rules", text: S.fmt(S.results.fixSynth, { reasons: plain.join("; ") }),
        buttons: [{ id: "newSeed", label: S.results.newSeed, newSeed: true }, { id: "keep", label: S.results.keep, keep: true }] };
    }
    var goal = Number(limitsUsed.cai_floor).toFixed(2);
    var out = { constraint: c.constraint, buttons: [] };
    var r = S.results;
    if (c.constraint === "window_energy") {
      var value = Math.floor(c.reached);
      out.text = c.blocked_by === "cai_floor"
        ? S.fmt(r.fixCai, { goal: goal, reached: Math.round(c.reached), limit: Math.round(c.limit) })
        : S.fmt(r.fixBudget, { reached: Math.round(c.reached), limit: Math.round(c.limit) });
      if (c.blocked_by === "cai_floor") out.buttons.push({ id: "lowerScore", label: r.lowerScore, override: { cai_floor: 0.85 } });
      out.buttons.push({ id: "allowFolds", label: S.fmt(r.allowFolds, { value: value }), override: { worst_window_energy: value } });
    } else if (c.constraint === "stem_length") {
      out.text = c.blocked_by === "cai_floor"
        ? S.fmt(r.fixStem, { goal: goal, reached: c.reached, limit: c.limit })
        : S.fmt(r.fixBudget, { reached: c.reached, limit: c.limit });
      if (c.blocked_by === "cai_floor") out.buttons.push({ id: "lowerScore", label: r.lowerScore, override: { cai_floor: 0.85 } });
      out.buttons.push({ id: "allowStem", label: S.fmt(r.allowStem, { value: c.reached }), override: { longest_allowed_stem: c.reached } });
    } else if (c.constraint === "start_region") {
      var sv = c.reached === null || c.reached === undefined ? null : Math.floor(c.reached);
      out.text = c.blocked_by === "cai_floor"
        ? S.fmt(r.fixStart, { goal: goal, reached: sv === null ? "not open" : sv, limit: Math.round(c.limit) })
        : S.fmt(r.fixBudget, { reached: sv === null ? "not open" : sv, limit: Math.round(c.limit) });
      if (c.blocked_by === "cai_floor") out.buttons.push({ id: "lowerScore", label: r.lowerScore, override: { cai_floor: 0.85 } });
      if (sv !== null) out.buttons.push({ id: "allowStart", label: S.fmt(r.allowStart, { value: sv }), override: { start_region_energy: sv } });
    } else if (c.constraint === "cai_floor") {
      out.text = S.fmt(r.fixScore, { reached: Number(c.reached).toFixed(2), limit: Number(c.limit).toFixed(2) });
      out.buttons.push({ id: "lowerScore", label: r.lowerScore, override: { cai_floor: 0.85 } });
    } else {
      out.text = S.fmt(r.fixOther, { reached: c.reached });
      if (c.blocked_by === "cai_floor") out.buttons.push({ id: "lowerScore", label: r.lowerScore, override: { cai_floor: 0.85 } });
    }
    out.buttons.push({ id: "newSeed", label: r.newSeed, newSeed: true });
    out.buttons.push({ id: "keep", label: r.keep, keep: true });
    return out;
  }

  /* DNA in blocks of 10 with running position numbers, 5 blocks per row. */
  function dnaRows(dna, perRow) {
    perRow = perRow || 5;
    var rows = [];
    for (var i = 0; i < dna.length; i += 10 * perRow) {
      var blocks = [];
      for (var j = i; j < Math.min(dna.length, i + 10 * perRow); j += 10) blocks.push(dna.substr(j, 10));
      rows.push({ start: i + 1, blocks: blocks });
    }
    return rows;
  }

  function fastaText(name, dna) {
    var lines = [">" + (name || "optimized_sequence")];
    for (var i = 0; i < dna.length; i += 60) lines.push(dna.substr(i, 60));
    return lines.join("\n") + "\n";
  }

  /* Flesch-Kincaid grade level over a list of strings. */
  function syllables(word) {
    word = word.toLowerCase().replace(/[^a-z]/g, "");
    if (!word) return 0;
    if (word.length <= 3) return 1;
    word = word.replace(/(?:[^laeiouy]es|ed|[^laeiouy]e)$/, "").replace(/^y/, "");
    var m = word.match(/[aeiouy]{1,2}/g);
    return Math.max(1, m ? m.length : 1);
  }
  function fleschKincaid(texts) {
    var sentences = 0, words = 0, syl = 0;
    texts.forEach(function (t) {
      var parts = String(t).split(/[.!?]+\s|[.!?]+$/).filter(function (p) { return p.trim(); });
      sentences += Math.max(1, parts.length);
      String(t).split(/\s+/).forEach(function (w) {
        var clean = w.replace(/[^A-Za-z']/g, "");
        if (clean) { words++; syl += syllables(clean); }
      });
    });
    if (!sentences || !words) return 0;
    return 0.39 * (words / sentences) + 11.8 * (syl / words) - 15.59;
  }

  var api = {
    AMINO: AMINO, DEFAULT_TEMPERATURE: DEFAULT_TEMPERATURE,
    analyzeProtein: analyzeProtein, removeInvalid: removeInvalid, translateFrame1: translateFrame1,
    validateTag: validateTag, validateSite: validateSite, validateStartDna: validateStartDna,
    validateTemperature: validateTemperature, validateSeed: validateSeed, validateRegions: validateRegions,
    detectRegions: detectRegions, hostTemperature: hostTemperature, buildRequests: buildRequests, mergeReports: mergeReports,
    statusCounts: statusCounts, bannerText: bannerText, checkResultText: checkResultText, techName: techName,
    fixPanel: fixPanel, dnaRows: dnaRows, fastaText: fastaText, fleschKincaid: fleschKincaid,
  };
  if (typeof module !== "undefined" && module.exports) module.exports = api;
  else root.PD_LOGIC = api;
})(typeof window !== "undefined" ? window : globalThis);
