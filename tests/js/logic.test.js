const test = require("node:test");
const assert = require("node:assert/strict");
const L = require("../../static/logic.js");
const S = require("../../static/strings.js");

const first = (a, id) => a.errors.concat(a.warnings, a.infos, a.notes).find((e) => e.id === id);

// ---- Errors section: one test per row, exact copy ---------------------------------
test("empty box", () => {
  const a = L.analyzeProtein("   ");
  assert.equal(a.blocked, true); assert.equal(a.empty, true);
  assert.equal(S.errors.empty, "Paste a protein sequence to start.");
});
test("letters that aren't amino acids", () => {
  const b = L.analyzeProtein("MKT-AY?K");
  assert.equal(first(b, "invalidLetters").text,
    "Found 2 letters that aren't amino acids (highlighted). Proteins use only A C D E F G H I K L M N P Q R S T V W Y.");
  assert.deepEqual(b.invalid.map((x) => x.char), ["-", "?"]);
  assert.equal(L.removeInvalid(b.sequence), "MKTAYK");
});
test("spaces, numbers, line breaks cleaned with one note", () => {
  const a = L.analyzeProtein("MKTAY 10\nIAKQR 20");
  assert.equal(a.sequence, "MKTAYIAKQR"); assert.equal(a.blocked, false);
  assert.equal(a.notes.filter((n) => n.id === "cleaned").length, 1);
  assert.equal(first(a, "cleaned").text, "We removed spaces and numbers from your sequence.");
});
test("FASTA header", () => {
  const a = L.analyzeProtein(">my_protein desc\nMKTAYIAKQRQISFVKSHFSRQ\nLEERLGLIEVQ");
  assert.equal(a.name, "my_protein desc");
  assert.equal(first(a, "fasta").text, "We used the name after > as your sequence name.");
  assert.equal(a.sequence, "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ");
  assert.equal(first(a, "cleaned"), undefined, "FASTA line wrapping is not 'cleaning'");
});
test("star at the end", () => {
  const a = L.analyzeProtein("MKTAYIAKQRQISFVKSHFSRQ*");
  assert.equal(a.sequence, "MKTAYIAKQRQISFVKSHFSRQ");
  assert.equal(first(a, "star").text, "We removed the * at the end. We add the stop signal for you.");
});
test("star in the middle blocks", () => {
  const a = L.analyzeProtein("MKTAYIAK*QRQISFVKSHFSRQ");
  assert.equal(a.blocked, true);
  assert.equal(first(a, "middleStar").text, "There's a * in the middle of your sequence. That would stop protein-making early.");
});
test("DNA pasted instead of protein", () => {
  const a = L.analyzeProtein("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG");
  assert.equal(a.dna, true); assert.equal(a.blocked, true);
  assert.equal(first(a, "dna").text, "This looks like DNA, not a protein. This tool starts from amino acids.");
  assert.equal(L.translateFrame1(a.dnaSequence), "MAIVMGR");
});
test("DNA detection needs >= 95% ACGTU", () => {
  assert.equal(L.analyzeProtein("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAGLLLL").dna, false); // 93%
  assert.equal(L.analyzeProtein("MKTAYIAKQRQISFVKSHFSRQ").dna, false);
});
test("X, B, Z, J, U, O", () => {
  const seq = "MKTAYIAKQRQISFVKSHFSRQ";
  assert.equal(first(L.analyzeProtein(seq + "X"), "badLetter").text,
    "X means an unknown amino acid, so we can't choose DNA for it. Replace it with a real one.");
  assert.equal(first(L.analyzeProtein(seq + "B"), "badLetter").text,
    "B means an ambiguous amino acid, so we can't choose DNA for it. Replace it with a real one.");
  for (const l of "ZJ") assert.match(first(L.analyzeProtein(seq + l), "badLetter").text, /ambiguous amino acid/);
  for (const l of "UO") assert.match(first(L.analyzeProtein(seq + l), "badLetter").text, /rare amino acid that we don't support/);
  assert.equal(L.analyzeProtein(seq + "X").blocked, true);
});
test("several FASTA records", () => {
  const a = L.analyzeProtein(">a\nMKTAYIAKQRQISFVKSHFSRQ\n>b\nMQIFVKTLTGKTITLEVEPSDT\n>c\nMKKKKKKKKKKKKKKKKKKKKKK");
  assert.equal(a.needsRecordChoice, true); assert.equal(a.blocked, true);
  assert.equal(first(a, "multi").text, "We found 3 sequences. Optimize them one at a time.");
  const b = L.analyzeProtein(">a\nMKTAYIAKQRQISFVKSHFSRQ\n>b\nMQIFVKTLTGKTITLEVEPSDT", { chosenRecord: 1 });
  assert.equal(b.sequence, "MQIFVKTLTGKTITLEVEPSDT"); assert.equal(b.blocked, false);
});
test("very short protein: warn under 20, block under 3", () => {
  assert.equal(first(L.analyzeProtein("MKTAYIAKQRQI"), "short").text,
    "That's only 12 amino acids. Most proteins are longer, so double-check it.");
  assert.equal(L.analyzeProtein("MKTAYIAKQRQI").blocked, false);
  assert.equal(L.analyzeProtein("MK").blocked, true);
});
test("too long", () => {
  const a = L.analyzeProtein("MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDG".repeat(80).slice(0, 3200), { maxLength: 1500 });
  assert.equal(a.blocked, true);
  assert.equal(first(a, "long").text, "This protein is 3,200 amino acids. The limit is 1,500. Try splitting it into domains.");
});
test("doesn't start with M is info only", () => {
  const a = L.analyzeProtein("KTAYIAKQRQISFVKSHFSRQ");
  assert.equal(a.blocked, false);
  assert.equal(first(a, "noMet").text,
    "Heads up: your protein doesn't start with M. We add a start signal (ATG) unless you tick 'My vector already provides the start codon'.");
});
test("no host picked", () => assert.equal(S.errors.noHost, "Pick at least one place to make your protein."));
test("bad tag letters", () => {
  assert.equal(L.validateTag("HHHHZH").error, "Tags use amino acid letters only. 'HHHHZH' has a Z.");
  assert.equal(L.validateTag("H".repeat(31)).ok, true);
  assert.match(L.validateTag("H".repeat(31)).warning, /Long tags/);
});
test("bad cut site", () => {
  assert.equal(L.validateSite("GATATX").error,
    "Cut sites use the letters A, C, G and T (plus standard ambiguity letters). 'GATATX' has an X.");
  assert.equal(L.validateSite("GAT").error, "Cut sites are usually 4 to 10 letters.");
  assert.equal(L.validateSite("GATATCGATATC").error, "Cut sites are usually 4 to 10 letters.");
  assert.equal(L.validateSite("gRatc").ok, true);
});
test("regions: range and overlap", () => {
  const e = L.validateRegions([{ start: 40, end: 52 }], 48);
  assert.equal(e[0], "This region ends at 52, but your protein is only 48 amino acids long.");
  const o = L.validateRegions([{ start: 1, end: 10 }, { start: 8, end: 20 }], 48);
  assert.deepEqual(o, ["These two regions overlap.", "These two regions overlap."]);
  assert.deepEqual(L.validateRegions([{ start: 1, end: 10 }, { start: 11, end: 20 }], 48), [null, null]);
});
test("start-of-gene DNA: letters, and U converted with a note", () => {
  assert.equal(L.validateStartDna("AGGAGX").error, "Use DNA letters only: A, C, G, T.");
  const r = L.validateStartDna("AGGAGGUAAA");
  assert.equal(r.dna, "AGGAGGTAAA"); assert.equal(r.note, "We changed U to T in your start-of-gene DNA.");
});
test("temperature range", () => {
  for (const bad of ["3.9", "45.1", "", "abc"]) assert.equal(L.validateTemperature(bad).error, "Choose a temperature between 4 and 45 °C.");
  assert.equal(L.validateTemperature("4").ok, true); assert.equal(L.validateTemperature("45").ok, true);
});
test("seed whole number", () => {
  assert.equal(L.validateSeed("4.5").error, "Use a whole number, like 42.");
  assert.equal(L.validateSeed("abc").error, "Use a whole number, like 42.");
  assert.equal(L.validateSeed("42").value, 42); assert.equal(L.validateSeed("").value, null);
});

// ---- Requests ---------------------------------------------------------------------
const DEFAULTS = { cai_floor: 0.9, worst_window_energy: -15, gc_min: 0.3, gc_max: 0.7, longest_allowed_stem: 7 };
const base = () => ({ protein: "MKTAYIAKQRQISFVKSHFSRQ", hosts: ["e_coli_bl21_de3"], goal: "BALANCED", seed: null,
  limits: {}, sites: [], regions: [], hedge: false, temperatureEdited: false });

test("default request sends only mode, protein, hosts", () => {
  const [r] = L.buildRequests(base(), DEFAULTS);
  assert.deepEqual(r, { protein: "MKTAYIAKQRQISFVKSHFSRQ", mode: "BALANCED", hosts: ["e_coli_bl21_de3"] });
});
test("only values that differ from defaults are sent", () => {
  const s = base(); s.limits = { cai: 0.9, mfe: -12, gcMin: 0.3, gcMax: 0.7, stem: 7 }; s.seed = 42; s.hedge = true;
  const [r] = L.buildRequests(s, DEFAULTS);
  assert.equal(r.worst_window_energy, -12);
  for (const k of ["cai_floor", "gc_min", "gc_max", "longest_allowed_stem"]) assert.equal(k in r, false, k);
  assert.equal(r.seed, 42); assert.equal(r.hedge, true);
});
test("M is added when missing", () => {
  const s = base(); s.protein = "KTAYIAKQRQISFVKSHFSRQ";
  assert.equal(L.buildRequests(s, DEFAULTS)[0].protein[0], "M");
});
test("suggested temperatures split hosts into requests; user-set temperature does not", () => {
  const s = base(); s.hosts = ["e_coli_bl21_de3", "s_cerevisiae", "human"];
  const reqs = L.buildRequests(s, DEFAULTS);
  assert.equal(reqs.length, 2);
  assert.deepEqual(reqs.map((r) => r.hosts).sort(), [["e_coli_bl21_de3", "human"], ["s_cerevisiae"]]);
  assert.equal(reqs.find((r) => r.hosts[0] === "s_cerevisiae").temperature, 30);
  assert.equal("temperature" in reqs.find((r) => r.hosts.length === 2), false);
  s.temperatureEdited = true; s.temperature = 25;
  const one = L.buildRequests(s, DEFAULTS); assert.equal(one.length, 1); assert.equal(one[0].temperature, 25);
});
test("regions only sent for folding-type goals; natural organism only when set", () => {
  const s = base(); s.regions = [{ start: 3, end: 9, type: "linker" }]; s.naturalOrganism = "human";
  let r = L.buildRequests(s, DEFAULTS)[0];
  assert.deepEqual(r.structural_regions, [{ start: 3, end: 9, kind: "linker" }]); assert.equal(r.native_organism, "human");
  s.goal = "PRODUCTION"; r = L.buildRequests(s, DEFAULTS)[0];
  assert.equal("structural_regions" in r, false); assert.equal("native_organism" in r, false);
});
test("vector ATG flag and tags", () => {
  const s = base(); s.vectorAtg = true; s.nTag = "HHHHHH"; s.cTag = "DYKDDDDK";
  const r = L.buildRequests(s, DEFAULTS)[0];
  assert.equal(r.vector_provides_start, true); assert.equal(r.n_tag, "HHHHHH"); assert.equal(r.c_tag, "DYKDDDDK");
});
test("tag presets are the spec's sequences", () => {
  const t = S.advanced.tags;
  assert.equal(t.his6.aa, "HHHHHH"); assert.equal(t.his6tev.aa, "HHHHHHENLYFQG");
  assert.equal(t.strep.aa, "WSHPQFEK"); assert.equal(t.flag.aa, "DYKDDDDK"); assert.equal(t.ha.aa, "YPYDVPDYA");
});
test("cut-site presets are the spec's sequences", () => {
  const e = S.advanced.enzymes;
  assert.deepEqual(S.advanced.presets.biobrick.map((n) => n + ":" + e[n]),
    ["EcoRI:GAATTC", "XbaI:TCTAGA", "SpeI:ACTAGT", "PstI:CTGCAG", "NotI:GCGGCCGC"]);
  assert.deepEqual(S.advanced.presets.common.map((n) => n + ":" + e[n]),
    ["NdeI:CATATG", "XhoI:CTCGAG", "BamHI:GGATCC", "NcoI:CCATGG", "HindIII:AAGCTT"]);
});
test("host table matches the spec (13 hosts, temperatures, 4 common)", () => {
  assert.equal(Object.keys(S.hosts).length, 13);
  assert.deepEqual(Object.keys(S.hosts).filter((h) => S.hosts[h].common), ["e_coli_bl21_de3", "e_coli_k12", "s_cerevisiae", "human"]);
  const temps = { e_coli_bl21_de3: 37, e_coli_k12: 37, s_cerevisiae: 30, human: 37, b_subtilis_168: 37, b_megaterium: 37,
    c_glutamicum: 30, l_lactis: 30, p_putida: 30, a_tumefaciens: 28, synechocystis_pcc6803: 30, p_pastoris: 30, c_reinhardtii: 25 };
  for (const [h, t] of Object.entries(temps)) assert.equal(S.hosts[h].temp, t, h);
});

// ---- Results ------------------------------------------------------------------------
const seq = (over) => Object.assign({ checks: [
  { id: "protein_identity", status: "pass", applicable: true, value: { matched: 359, total: 359 } },
  { id: "hairpins", status: "fail", applicable: true, value: { longest_stem: 9, strongest: -19.2 } },
  { id: "terminators", status: "pass", applicable: false, value: { found: 0 } } ], conflicts: [] }, over);
test("banner wording and counts ignore non-applicable checks", () => {
  const b = L.bannerText(seq());
  assert.equal(b.ok, false); assert.equal(b.text, "Your DNA is ready, but 1 check needs a look.");
  const ok = L.bannerText({ checks: seq().checks.map((c) => Object.assign({}, c, { status: "pass" })) });
  assert.equal(ok.text, "Your DNA is ready. All 2 checks passed.");
});
test("check result sentences", () => {
  assert.equal(L.checkResultText({ id: "protein_identity", status: "pass", value: { matched: 359, total: 359 } }), "359 of 359 amino acids match");
  assert.equal(L.checkResultText({ id: "rare_codons", status: "pass", value: { found: 0, lowest: 0.339, cutoff: 0.3 } }), "0 found; lowest score 0.34");
  assert.equal(L.checkResultText({ id: "codon_score", status: "pass", value: { cai: 0.931, goal: 0.9 } }), "0.93 (goal 0.90 or higher)");
  assert.equal(L.checkResultText({ id: "hairpins", status: "pass", value: { longest_stem: 6, strongest: -12.4 } }), "Longest 6 letters; strongest -12.4");
  // a value just past the -15 limit must not display as if it were exactly on it
  assert.equal(L.checkResultText({ id: "hairpins", status: "fail", value: { longest_stem: 7, strongest: -15.4 } }), "Longest 7 letters; strongest -15.4");
  assert.equal(L.checkResultText({ id: "synthesis", status: "fail", value: { gc_min: 0.296, gc_max: 0.6 } }), "Repeats or lopsided DNA; GC 29.6-60.0%");
  assert.equal(L.checkResultText({ id: "start_region", status: "review", value: {} }).startsWith("Not checked"), true);
  assert.equal(L.techName({ id: "rare_codons", value: { cutoff: 0.3 } }), "Rare codons (score under 0.3)");
});
test("fix panel names the numbers and only offers explicit fixes", () => {
  const p = L.fixPanel(seq({ conflicts: [{ constraint: "window_energy", limit: -15, reached: -18.6, blocked_by: "cai_floor" }] }), { cai_floor: 0.9 });
  assert.equal(p.text, "To reach a codon match score of 0.90, the DNA has to use choices that make the gene's message fold more than your limit (strongest fold -19, limit -15).");
  assert.deepEqual(p.buttons.map((b) => b.label), ["Lower the score goal to 0.85 and re-run", "Allow stronger folds (-19) and re-run", "Keep this result"]);
  assert.deepEqual(p.buttons[0].override, { cai_floor: 0.85 }); assert.deepEqual(p.buttons[1].override, { worst_window_energy: -19 });
  assert.equal(L.fixPanel(seq({ conflicts: [] }), { cai_floor: 0.9 }), null);
  assert.equal(L.fixPanel(seq({ conflicts: [{ constraint: "sequence_rules", reasons: ["x"] }] }), { cai_floor: 0.9 }), null);
});
test("DNA blocks of 10 with positions; FASTA wraps at 60", () => {
  const dna = "ACGT".repeat(30);
  const rows = L.dnaRows(dna); assert.equal(rows[0].blocks[0].length, 10); assert.equal(rows[1].start, 51);
  assert.equal(L.fastaText("x", dna).split("\n")[1].length, 60);
});

// ---- Acceptance: reading level, glossary, performance ---------------------------------
test("default-view copy reads at grade 9 or lower (Flesch-Kincaid)", () => {
  const grade = L.fleschKincaid(S.defaultViewText());
  assert.ok(grade <= 9, "grade " + grade.toFixed(1));
});
test("no raw internal IDs or JSON on the default view text", () => {
  const all = S.defaultViewText().join(" ");
  for (const id of Object.keys(S.hosts)) assert.equal(all.includes(id), false, id);
  assert.equal(/[{}\[\]]/.test(all.replace(/\{n\}/g, "")), false);
});
test("glossary defines every term the checklist and default view rely on", () => {
  const keys = S.glossary.terms.map((t) => t[0]);
  for (const k of ["host", "codon", "cai", "rare", "hairpin", "fold", "rbs", "atg", "stop", "cutsite", "goldengate", "biobrick",
    "tag", "tev", "linker", "folding", "denovo", "gc", "repeat", "synthesis", "backup", "seed", "fasta"]) assert.ok(keys.includes(k), k);
  assert.equal(S.glossary.terms.length, 23);
  for (const [, , text] of S.glossary.terms) assert.ok(text.split(/[.?!]\s/).length <= 2, "tooltip too long: " + text);
});
test("validating 2,000 amino acids takes under 100 ms", () => {
  const p = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQ".repeat(40).slice(0, 2000);
  L.analyzeProtein(p); // warm up
  const t = process.hrtime.bigint(); L.analyzeProtein(p);
  const ms = Number(process.hrtime.bigint() - t) / 1e6;
  assert.ok(ms < 100, ms + " ms");
});
test("auto-detected linkers", () => {
  assert.deepEqual(L.detectRegions("MKAAAGGGGSGGGGSAAAK"), [{ start: 6, end: 15, type: "linker", auto: true }]);
});
