/* All user-visible text lives here so teams can translate it.
   Nothing else in the client contains display copy. Pure data: usable in the
   browser (window.PD_STRINGS) and in node tests (require). */
(function (root) {
  "use strict";

  var S = {};

  // Text that appears on the default view (used by the reading-level check).
  S.app = {
    title: "Codon Optimizer",
    intro:
      "Turn a protein into DNA that cells can read well. Paste your protein, pick where it will be made, and we'll check the result.",
    tryExample: "Try an example",
    clear: "Clear",
    optimize: "Optimize my sequence",
    optimizeDisabledTip: "Paste a protein sequence to start.",
    belowButton: "We check your DNA for common problems before you order it.",
    glossaryLink: "Glossary",
    advancedLink: "Advanced settings",
    advancedChanged: "Advanced settings ({n} changed)",
    privacy: "We don't save your sequence. It is processed on our server and not stored.",
    exampleName: "Ubiquitin (example)",
    // A short, familiar, harmless protein.
    exampleProtein:
      "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
  };

  S.steps = {
    protein: {
      label: "1. Paste your protein sequence.",
      help: "Use one-letter amino acid codes, like MKTAY. You can paste a FASTA file too.",
      counter: "{n} amino acids",
      counterOne: "1 amino acid",
    },
    host: {
      label: "2. Where will your protein be made?",
      help: "This is the living cell that reads your DNA. Not sure? Keep E. coli.",
      more: "More organisms",
      search: "Search organisms",
      searchHelp: "Try yeast, coli or Bacillus.",
      noMatch: "No organism matches that search.",
      tabsNote: "Results will appear in tabs.",
    },
    goal: {
      label: "3. What matters most?",
    },
    natural: {
      label: "Is this protein found in nature?",
      no: "No, I designed it",
      yes: "Yes",
      organism: "Which organism is it from?",
      organismHelp: "Pick the closest match from the list.",
    },
  };

  S.goals = {
    PRODUCTION: { title: "Make the most protein", desc: "Best when you need a lot of protein, for example to purify it.", tag: "PRODUCTION" },
    FOLDING: { title: "Make protein that folds correctly", desc: "Best for delicate proteins that must fold just right to work.", tag: "FOLDING" },
    BALANCED: { title: "Best balance", desc: "A good mix of both. Recommended for most projects.", tag: "BALANCED", recommended: "Recommended" },
    ALL: { title: "Show me all three", desc: "See each version side by side.", tag: "ALL" },
  };

  // Display names for result titles (mode as returned by the engine).
  S.goalNames = {
    PRODUCTION: "Most protein",
    FOLDING: "Folds correctly",
    BALANCED: "Best balance",
    "BALANCED (diversified)": "Backup version",
  };

  S.hosts = {
    e_coli_bl21_de3: { short: "E. coli (BL21)", title: "E. coli, protein-making strain (BL21)", latin: "Escherichia coli BL21", desc: "Bacteria. The workhorse strain for making lots of protein.", common: true, temp: 37 },
    e_coli_k12: { short: "E. coli (K-12)", title: "E. coli, standard lab strain (K-12)", latin: "Escherichia coli K-12", desc: "Bacteria. The everyday lab strain, often used for cloning.", common: true, temp: 37 },
    s_cerevisiae: { title: "Baker's yeast", latin: "Saccharomyces cerevisiae", desc: "Yeast. A familiar single-celled fungus and popular lab host.", common: true, temp: 30 },
    human: { title: "Human cells", latin: "Homo sapiens", desc: "Human cell lines, for example HEK293.", common: true, temp: 37 },
    b_subtilis_168: { title: "Bacillus subtilis 168", latin: "Bacillus subtilis 168", desc: "Bacteria. A soil bacterium that can release proteins outside the cell.", temp: 37 },
    b_megaterium: { title: "Bacillus megaterium", latin: "Bacillus megaterium", desc: "Bacteria. A large soil bacterium used in industry to make proteins.", temp: 37 },
    c_glutamicum: { title: "Corynebacterium glutamicum", latin: "Corynebacterium glutamicum", desc: "Bacteria. Used in industry to make amino acids.", temp: 30 },
    l_lactis: { title: "Lactococcus lactis", latin: "Lactococcus lactis", desc: "Bacteria. A dairy bacterium, common in food uses.", temp: 30 },
    p_putida: { title: "Pseudomonas putida", latin: "Pseudomonas putida", desc: "Bacteria. A tough soil bacterium that breaks down many chemicals.", temp: 30 },
    a_tumefaciens: { title: "Agrobacterium tumefaciens", latin: "Agrobacterium tumefaciens", desc: "Bacteria. Naturally moves DNA into plants, used to engineer them.", temp: 28 },
    synechocystis_pcc6803: { title: "Synechocystis PCC 6803", latin: "Synechocystis sp. PCC 6803", desc: "Cyanobacteria. A photosynthetic microbe that runs on light.", temp: 30 },
    p_pastoris: { title: "Pichia pastoris", latin: "Komagataella phaffii", desc: "Yeast. Used to make and release large amounts of protein.", temp: 30 },
    c_reinhardtii: { title: "Chlamydomonas reinhardtii", latin: "Chlamydomonas reinhardtii", desc: "Algae. A green alga that grows on light.", temp: 25 },
  };
  S.defaultHost = "e_coli_bl21_de3";

  // Technical labels for the collapsed "Full technical report".
  S.tech = {
    cai: "CAI", tai: "tAI", taiNone: "Not available for this host", minW: "Lowest codon weight (w)",
    belowW: "Codons under the rare cutoff", startDg: "Start-region fold energy (-20 to +40)",
    startOpen: "Start region single-stranded (-15 to +20)", worstMfe: "Worst 60-letter window fold energy",
    stem: "Longest stem (contiguous pairs)", terminators: "Terminator-like motifs", inverted: "Inverted repeats (8+)",
    gc: "GC overall", gcWindow: "GC, 50-letter windows (lowest to highest)", repeats: "Repeated 15-letter pieces",
    pauses: "Slow spots placed", yes: "Yes", no: "No", none: "Not available", kcal: "kcal/mol",
    sources: "Codon table source", notes: "Notes from the engine",
  };

  S.advanced = {
    reset: "Reset to defaults",
    developer: "Developer options",
    groups: {
      vector: "Your vector",
      tags: "Tags",
      sites: "Cut sites",
      regions: "Protein regions",
      growth: "Growth",
      results: "Results",
      limits: "Expert limits",
    },
    start: {
      label: "Start of your gene",
      help: "The short DNA just before your gene that tells the cell where to start reading. We use it to check the start of your gene for problems.",
      skip: "Not sure / skip",
      paste: "Paste my own DNA",
      pasteLabel: "Your start-of-gene DNA",
      partsGroup: "Registry RBS parts for your host",
      partsNone: "No Registry RBS parts are listed for this host yet.",
    },
    vectorAtg: {
      label: "My vector already provides the start codon (ATG)",
      help: "Tick this if your plasmid adds the ATG for you.",
    },
    nTag: { label: "Add a tag to the front (N-terminal)", help: "A short extra piece on your protein that helps you purify or detect it." },
    cTag: { label: "Add a tag to the end (C-terminal)", help: "If you add an end tag, the stop signal goes after it." },
    tagCustom: "Custom",
    tagCustomLabel: "Custom tag (amino acids)",
    tagNone: "None",
    tags: {
      none: { name: "None", aa: "" },
      his6: { name: "6xHis", aa: "HHHHHH" },
      his6tev: { name: "6xHis + TEV cut site", aa: "HHHHHHENLYFQG" },
      strep: { name: "Strep-tag II", aa: "WSHPQFEK" },
      flag: { name: "FLAG", aa: "DYKDDDDK" },
      ha: { name: "HA", aa: "YPYDVPDYA" },
      custom: { name: "Custom", aa: "" },
    },
    sites: {
      label: "DNA cut sites to avoid",
      help: "We keep these out of your gene so your cloning enzymes don't chop it up.",
      locked: "Always avoided: BsaI, BsmBI, BbsI, SapI",
      biobrick: "iGEM BioBrick (EcoRI, XbaI, SpeI, PstI, NotI)",
      common: "Common cloning (NdeI, XhoI, BamHI, NcoI, HindIII)",
      addOwn: "Add my own cut site",
      addName: "Enzyme name",
      addSite: "DNA site",
      add: "Add",
      remove: "Remove {name}",
      suggestions: "Suggestions",
    },
    // Real recognition sequences for the two presets and the search list.
    enzymes: {
      EcoRI: "GAATTC", XbaI: "TCTAGA", SpeI: "ACTAGT", PstI: "CTGCAG", NotI: "GCGGCCGC",
      NdeI: "CATATG", XhoI: "CTCGAG", BamHI: "GGATCC", NcoI: "CCATGG", HindIII: "AAGCTT",
      EcoRV: "GATATC", KpnI: "GGTACC", SacI: "GAGCTC", SalI: "GTCGAC", SphI: "GCATGC", ClaI: "ATCGAT",
    },
    presets: {
      biobrick: ["EcoRI", "XbaI", "SpeI", "PstI", "NotI"],
      common: ["NdeI", "XhoI", "BamHI", "NcoI", "HindIII"],
    },
    regions: {
      label: "Flexible parts of your protein",
      help: "Places where slowing the cell down a little can help your protein fold. We find likely ones for you.",
      start: "Start", end: "End", type: "Type",
      linker: "Linker", loop: "Loop", domain_boundary: "Domain boundary",
      add: "Add region", remove: "Remove", reset: "Reset to auto-detected",
      none: "No flexible parts found. You can add some.",
      highlighted: "Your sequence with likely linkers highlighted",
      relevant: "Only used for the folding, balance and all-three goals.",
    },
    temperature: {
      label: "Growth temperature (°C)",
      help: "The temperature your cells grow at. We use it to predict whether your gene's message folds up in a way that could block protein production.",
      quick: "Quick choices",
      hostNote: "Suggested for your host. Please confirm it.",
      mixed: "Your hosts grow at different temperatures, so each one uses its own suggested value. Type a number to use one temperature for all.",
    },
    backup: { label: "Also give me a backup version", help: "We'll make a second version with very different DNA choices, in case the first one doesn't work well." },
    seed: { label: "Repeatable results", help: "Enter a number to get the same result every time. Leave blank for a fresh one.", used: "Number used: {seed}. Enter it to get this exact result again." },
    limits: {
      warning: "Most teams should leave these alone. Changing them can make some checks fail.",
      cai: "Codon match score goal",
      mfe: "Worst allowed RNA fold (kcal/mol)",
      gcMin: "Lowest GC (%)",
      gcMax: "Highest GC (%)",
      stem: "Longest allowed stem (letters)",
    },
    developerOptions: {
      export: "Export regions as JSON",
      import: "Import regions from JSON",
      importLabel: "Regions JSON",
      importApply: "Use these regions",
      copyRequest: "Copy API request",
      copied: "Copied.",
      importBad: "That JSON isn't a list of {start, end, kind} items.",
    },
  };

  // Exact-copy messages. {placeholders} are filled by the client.
  S.errors = {
    empty: "Paste a protein sequence to start.",
    invalidLetters: "Found {n} letters that aren't amino acids (highlighted). Proteins use only A C D E F G H I K L M N P Q R S T V W Y.",
    invalidLettersOne: "Found 1 letter that isn't an amino acid (highlighted). Proteins use only A C D E F G H I K L M N P Q R S T V W Y.",
    removeThem: "Remove them",
    cleaned: "We removed spaces and numbers from your sequence.",
    fastaHeader: "We used the name after > as your sequence name.",
    trailingStar: "We removed the * at the end. We add the stop signal for you.",
    middleStar: "There's a * in the middle of your sequence. That would stop protein-making early.",
    dna: "This looks like DNA, not a protein. This tool starts from amino acids.",
    translate: "Translate it for me",
    badLetterX: "{letter} means {kind}, so we can't choose DNA for it. Replace it with a real one.",
    kindX: "an unknown amino acid",
    kindAmbiguous: "an ambiguous amino acid",
    kindRare: "a rare amino acid that we don't support",
    multi: "We found {n} sequences. Optimize them one at a time.",
    useFirst: "Use the first one",
    chooseOne: "Choose which one",
    chooseLabel: "Sequence to use",
    tooShortWarn: "That's only {n} amino acids. Most proteins are longer, so double-check it.",
    tooShortBlock: "That's only {n} amino acids. We need at least 3 to make DNA.",
    tooLong: "This protein is {n} amino acids. The limit is {limit}. Try splitting it into domains.",
    noMet: "Heads up: your protein doesn't start with M. We add a start signal (ATG) unless you tick 'My vector already provides the start codon'.",
    noHost: "Pick at least one place to make your protein.",
    organismRequired: "Choose the organism this protein comes from, or pick 'No, I designed it'.",
    fixFields: "Fix the problems shown, then try again.",
    tagBad: "Tags use amino acid letters only. '{tag}' has a {letter}.",
    tagLong: "That tag is {n} amino acids. Long tags can change how your protein folds.",
    siteBad: "Cut sites use the letters A, C, G and T (plus standard ambiguity letters). '{site}' has an {letter}.",
    siteLength: "Cut sites are usually 4 to 10 letters.",
    siteName: "Give this cut site a name.",
    regionRange: "This region ends at {end}, but your protein is only {n} amino acids long.",
    regionOrder: "The start of a region must come before its end.",
    regionStart: "Regions start at position 1 or higher.",
    regionOverlap: "These two regions overlap.",
    startDna: "Use DNA letters only: A, C, G, T.",
    startDnaU: "We changed U to T in your start-of-gene DNA.",
    temperature: "Choose a temperature between 4 and 45 °C.",
    seed: "Use a whole number, like 42.",
    number: "Enter a number.",
    server: "Something went wrong on our side. Your settings are saved, so you can try again.",
    offline: "You're offline. Reconnect and try again.",
    tryAgain: "Try again",
    removeRegion: "Remove region {n}",
    copyError: "Copy error details",
    problems: "{n} problems found",
    problemsOne: "1 problem found",
  };

  S.loading = {
    title: "Working on it",
    steps: ["Choosing DNA for each amino acid", "Removing rare choices", "Checking for hairpins", "Running final checks"],
    elapsed: "{s} seconds so far",
    cancel: "Cancel",
    slow: "This is taking longer than usual. Long proteins can take a few minutes. Keep waiting or cancel.",
    cancelled: "Cancelled. Your settings are saved.",
  };

  S.results = {
    ready: "Your DNA is ready. All {n} checks passed.",
    readySome: "Your DNA is ready, but {k} checks need a look.",
    readyOne: "Your DNA is ready, but 1 check needs a look.",
    statusPass: "Pass",
    statusReview: "Take a look",
    statusFail: "Needs fixing",
    copyDna: "Copy DNA",
    copied: "Copied.",
    downloadFasta: "Download FASTA",
    showProtein: "Show the protein this makes",
    hideProtein: "Hide the protein",
    proteinCaption: "The protein your DNA makes",
    dnaCaption: "Your DNA, in blocks of 10 letters",
    summary: "Summary",
    summaryHost: "Host",
    summaryGoal: "Goal",
    summaryChecks: "Checks passed",
    ofChecks: "{p} of {t}",
    mustPass: "Must pass",
    quality: "Quality",
    whatThisMeans: "What this means",
    editSettings: "Edit settings",
    runAgain: "Run again with the same seed",
    fullReport: "Full technical report",
    startNote: "Your gene starts without the ATG because your vector adds it.",
    regionsNote: "We used {n} flexible parts we found in your protein.",
    fixTitle: "We couldn't meet all your goals at once.",
    fixCai: "To reach a codon match score of {goal}, the DNA has to use choices that make the gene's message fold more than your limit (strongest fold {reached}, limit {limit}).",
    fixStem: "To reach a codon match score of {goal}, the DNA has to use choices that make a hairpin longer than your limit (longest {reached} letters, limit {limit}).",
    fixStart: "To reach a codon match score of {goal}, the start of your gene folds up more than your limit (start fold {reached}, limit {limit}).",
    fixBudget: "We tried several times but couldn't fully meet this limit (reached {reached}, limit {limit}).",
    fixScore: "The DNA could only reach a codon match score of {reached}, and your goal is {limit}.",
    fixSynth: "We couldn't remove every problem that makes this DNA hard to make: {reasons}. Another random try often clears it.",
    synthReasons: { "GC window high": "a stretch with too much G and C", "GC window low": "a stretch with too little G and C",
      "repeated 15-mer": "a piece that repeats", "homopolymer run": "a long run of one letter" },
    synthOther: "a cut site or other rule",
    newSeed: "Try again with a new random seed",
    fixOther: "We couldn't remove every problem in this check ({reached} left).",
    lowerScore: "Lower the score goal to 0.85 and re-run",
    allowFolds: "Allow stronger folds ({value}) and re-run",
    allowStem: "Allow longer hairpins ({value} letters) and re-run",
    allowStart: "Allow a less open start ({value}) and re-run",
    keep: "Keep this result",
    kept: "Kept. This result still has the problem described above.",
    copyFailed: "We couldn't copy that automatically. Select the text and copy it yourself.",
  };

  // Checklist rows: plain name, technical name, result templates, meaning.
  S.parts = {
    title: "Find parts for this gene",
    intro: "We look up promoters, ribosome binding sites and terminators from the iGEM Registry that pass the BioBrick (RFC10) check.",
    honest: "We can't predict how much protein a part makes. The Registry has no measured strength for almost all parts, so we don't guess.",
    levelLabel: "How much protein do you want?",
    levels: { any: "No preference", low: "Less", moderate: "A medium amount", high: "More" },
    levelHelp: "Only used when a part's own description says strong, weak or medium.",
    find: "Find parts",
    finding: "Looking up parts",
    kinds: { promoter: "Promoters", rbs: "Ribosome binding sites", terminator: "Terminators" },
    counts: "{p} of {c} parts passed our checks",
    none: "No parts passed the checks for this host.",
    available: "In the Registry stock",
    notAvailable: "Not currently in stock",
    hostTagged: "Made for your host",
    hostUnknown: "Host not recorded",
    why: "Why this part",
    viewRegistry: "View in the Registry",
    strengthNone: "No measured strength on file",
    regulation: { constitutive: "Always on (no inducer needed)", regulated: "Needs an inducer or repressor to switch on", unknown: "Switching not recorded" },
    strengthKnown: "Measured strength: {value} {unit}",
    hintText: "Description says {word} (not a measurement)",
    cdsBad: "Your gene has a cut site that RFC10 doesn't allow ({sites}). Fix it before assembling.",
    cdsOk: "Your gene passes the RFC10 check.",
    unsupported: "We don't have Registry host tags for this organism yet, so we can't match parts to it.",
    limits: "As soon as parts are found, we join the top-ranked ones with your gene below and check the whole plasmid, including the joins.",
    source: "Registry data: {file}, checksum {sha}",
    error: "Couldn't look up parts. Try again in a moment.",
  };

  S.build = {
    title: "Build the full plasmid",
    intro: "Choose one promoter, one ribosome binding site and one terminator above. We join them with your gene and check the whole plasmid.",
    use: "Use this part",
    button: "Build the plasmid",
    building: "Building",
    needParts: "Find parts first, then pick one of each kind.",
    backbone: "Backbone: {name} ({bp} letters). The Registry describes it as: {desc}.",
    junctionOk: "No illegal cut sites in the whole plasmid, including where the parts join.",
    junctionBad: "The finished plasmid has a cut site the BioBrick standard doesn't allow.",
    violation: "{enzyme} site at positions {start}-{end}{where}",
    atJoin: " (right at a join between parts)",
    flankNote: "The backbone's own EcoRI, XbaI, SpeI, PstI and NotI sites at its ends are expected and ignored.",
    sizeTitle: "Size",
    size: "{total} letters in total: {backbone} backbone plus {insert} insert.",
    sizeWithin: "Within the size range studied.",
    sizeBeyond: "Larger than the size range studied.",
    gcTitle: "GC balance of the insert",
    gc: "{gc}% GC overall. Its most extreme 50-letter windows are {lo}% and {hi}%.",
    gcTwistBad: "A 50-letter window is below 10% or above 90% GC. Twist calls that high-complexity, which is harder to make.",
    gcTwistOk: "No 50-letter window is below 10% or above 90% GC (the limits Twist states).",
    gcProject: "Overall GC is outside 30-70%. That range is this project's own guide, not a company limit.",
    synthBad: "The insert is longer than {limit} letters, the longest single gene Twist makes.",
    layoutTitle: "What is in the plasmid",
    scarNote: "Scars: {std} between most parts and {rbs} between the ribosome binding site and your gene.",
    seqTitle: "Your finished plasmid",
    copy: "Copy plasmid DNA",
    download: "Download FASTA",
    copied: "Copied.",
    sources: "Sources",
    sourceLabels: ["Registry: RFC10 standard", "Registry: scars", "Twist: sequence limits"],
    error: "Couldn't build the plasmid. Try again in a moment.",
    kind: { promoter: "Promoter", rbs: "Ribosome binding site", terminator: "Terminator", cds: "Your gene", scar: "Scar" },
    readyTitle: "Ready to order?",
    readyYes: "Yes — nothing here blocks ordering this construct as designed.",
    readyNo: "Not yet — this construct has a problem to fix first:",
    autoNote: "Built automatically from the top-ranked part in each list above. Pick different parts and press “Build the plasmid” to re-check.",
  };

  S.checks = {
    protein_identity: {
      name: "Your protein is made correctly", tech: "Protein identity",
      pass: "{matched} of {total} amino acids match", fail: "Only {matched} of {total} amino acids match",
      means: "If we read this DNA back, we get exactly your protein.",
    },
    start_stop: {
      name: "Start and stop signals", tech: "Stop codon, no early stops",
      pass: "Ends with a stop; no early stops", fail: "Missing a stop or stops too early",
      means: "The cell knows where to stop and won't stop too soon.",
    },
    cut_sites: {
      name: "No accidental cut sites", tech: "Restriction sites",
      pass: "0 found", fail: "{found} found",
      means: "Your cloning enzymes won't chop up your gene.",
    },
    rare_codons: {
      name: "No rare DNA words", tech: "Rare codons (score under {cutoff})",
      pass: "0 found; lowest score {lowest}", fail: "{found} found; lowest score {lowest}",
      means: "Rare choices can make the cell stall while building your protein.",
    },
    codon_score: {
      name: "Codon match score", tech: "CAI",
      pass: "{cai} (goal {goal} or higher)", fail: "{cai} (goal {goal} or higher)",
      means: "How closely the DNA uses the host's favorite choices. Higher usually means more protein.",
    },
    synthesis: {
      name: "Easy for companies to make", tech: "Repeats, long runs, GC balance",
      pass: "No repeats; GC {gcMin}-{gcMax}%", fail: "Repeats or lopsided DNA; GC {gcMin}-{gcMax}%",
      means: "Repeated or lopsided DNA is harder, and sometimes impossible, to synthesize.",
    },
    start_region: {
      name: "Gene start is open", tech: "Start-region folding",
      pass: "Open ({dG} kcal/mol)", fail: "Folded up ({dG} kcal/mol)",
      review: "Not checked. Add your start-of-gene DNA in Advanced settings.",
      means: "If the start of the message folds up, the cell struggles to begin.",
    },
    hairpins: {
      name: "No tight hairpins", tech: "Stem length and window energy",
      pass: "Longest {stem} letters; strongest {strongest}", fail: "Longest {stem} letters; strongest {strongest}",
      means: "Hairpins are spots where the message sticks to itself and can slow production.",
    },
    terminators: {
      name: "No stop-reading signals (bacteria)", tech: "Terminator-like motifs",
      pass: "None", fail: "{found} found",
      means: "Some patterns tell bacteria to stop reading early.",
    },
    mirror_repeats: {
      name: "No mirror repeats", tech: "Inverted repeats",
      pass: "None", fail: "{found} found",
      means: "Mirror-image stretches can fold into hairpins and confuse synthesis.",
    },
  };

  S.glossary = {
    title: "Glossary",
    intro: "Plain-English meanings for the terms used in this tool.",
    back: "Back to the optimizer",
    help: "What does {term} mean?",
    close: "Close",
    terms: [
      ["host", "Host", "The living cell that reads your DNA and builds your protein, such as E. coli or yeast."],
      ["codon", "Codon", "A three-letter DNA “word” that tells the cell which amino acid to add next. Most amino acids can be written several ways."],
      ["cai", "Codon match score (CAI)", "How closely your DNA uses the host's favorite word choices, from 0 to 1. Higher usually means more protein."],
      ["rare", "Rare codon", "A word choice the host hardly ever uses. Too many can make the cell slow down or stall."],
      ["hairpin", "Hairpin", "A spot where the gene's message folds back and sticks to itself. Tight hairpins can slow or block protein-making."],
      ["fold", "Fold energy (kcal/mol)", "How tightly the message sticks to itself. More negative means stickier, so closer to zero is better here."],
      ["rbs", "Start of your gene (RBS or Kozak)", "The short DNA just before your gene that tells the cell where to start reading."],
      ["atg", "Start codon (ATG)", "The three letters that mark where protein-building begins."],
      ["stop", "Stop codon", "Three letters that tell the cell to stop building the protein."],
      ["cutsite", "Cut site (restriction site)", "A short DNA pattern that a cloning enzyme cuts. We keep them out of your gene so it isn't chopped up."],
      ["goldengate", "Golden Gate", "A cloning method that uses enzymes such as BsaI and BsmBI to join DNA pieces in order."],
      ["biobrick", "BioBrick (RFC[10])", "The classic iGEM part standard. It uses EcoRI, XbaI, SpeI and PstI cut sites for joining parts."],
      ["tag", "Tag", "A short extra piece added to your protein so you can purify or detect it, such as a 6xHis tag."],
      ["tev", "TEV cut site", "A short amino acid sequence that lets you snip a tag off later."],
      ["linker", "Linker", "A flexible stretch of amino acids, often glycine and serine, that joins two parts of a protein."],
      ["folding", "Folding", "How a protein twists into its working 3D shape."],
      ["denovo", "De novo protein", "A protein you designed from scratch, not one found in nature."],
      ["gc", "GC content", "The share of DNA letters that are G or C. Very high or very low values make DNA harder to make."],
      ["repeat", "Repeat", "The same stretch of DNA appearing more than once. Repeats make DNA harder to synthesize and clone."],
      ["synthesis", "Synthesis-friendly", "DNA that a company can make reliably, without long repeats, extreme GC or long runs of one letter."],
      ["backup", "Backup version", "A second version with very different DNA choices, in case the first one doesn't work well."],
      ["seed", "Repeatable results (seed)", "A number that makes the tool give the same answer each time. Leave blank for a new one."],
      ["fasta", "FASTA", "A simple text format for sequences. A name line starting with > is followed by the letters."],
    ],
  };

  // Keys whose text is on the default view. Used by the reading-level check.
  S.defaultViewText = function () {
    var out = [];
    var push = function (x) { if (typeof x === "string" && x) out.push(x); };
    push(S.app.intro); push(S.app.belowButton); push(S.app.privacy); push(S.app.optimizeDisabledTip);
    push(S.steps.protein.label); push(S.steps.protein.help);
    push(S.steps.host.label); push(S.steps.host.help);
    push(S.steps.goal.label); push(S.steps.natural.label);
    ["PRODUCTION", "FOLDING", "BALANCED", "ALL"].forEach(function (g) { push(S.goals[g].title); push(S.goals[g].desc); });
    Object.keys(S.hosts).forEach(function (h) { if (S.hosts[h].common) { push(S.hosts[h].title); push(S.hosts[h].desc); } });
    return out;
  };

  S.fmt = function (template, vars) {
    return String(template).replace(/\{(\w+)\}/g, function (m, k) {
      return vars && vars[k] !== undefined && vars[k] !== null ? String(vars[k]) : m;
    });
  };

  if (typeof module !== "undefined" && module.exports) module.exports = S;
  else root.PD_STRINGS = S;
})(typeof window !== "undefined" ? window : globalThis);
