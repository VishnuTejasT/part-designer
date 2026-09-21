/* Codon Optimizer client. Builds the form once, then updates pieces in place so
   typing never loses focus. All text comes from strings.js; all rules from logic.js. */
(function () {
  "use strict";
  var S = window.PD_STRINGS, L = window.PD_LOGIC, F = S.fmt;
  var HOST_IDS = Object.keys(S.hosts);

  /* ---------- tiny DOM helpers ---------- */
  function h(tag, attrs) {
    var el = document.createElement(tag);
    Object.keys(attrs || {}).forEach(function (k) {
      var v = attrs[k];
      if (v === false || v === null || v === undefined) return;
      if (k === "class") el.className = v;
      else if (k === "text") el.textContent = v;
      else if (k.slice(0, 2) === "on") el.addEventListener(k.slice(2), v);
      else if (k === "value") el.value = v;
      else el.setAttribute(k, v === true ? "" : v);
    });
    for (var i = 2; i < arguments.length; i++) add(el, arguments[i]);
    return el;
  }
  function add(el, c) {
    if (c === null || c === undefined || c === false) return;
    if (Array.isArray(c)) c.forEach(function (x) { add(el, x); });
    else el.appendChild(typeof c === "string" ? document.createTextNode(c) : c);
  }
  function clear(el) { while (el.firstChild) el.removeChild(el.firstChild); return el; }
  var uid = 0;
  function id(prefix) { uid += 1; return prefix + uid; }
  function hostName(hid) { var x = S.hosts[hid]; return x.short || x.title; }

  /* ---------- state ---------- */
  var DEFAULT_LIMITS = { max_protein_length: 1500, temperature_range: [4, 45],
    defaults: { temperature: 37, cai_floor: 0.9, rare_codon_cutoff: 0.3, gc_min: 0.3, gc_max: 0.7,
      longest_allowed_stem: 7, worst_window_energy: -15, start_region_energy: -5 } };
  var api = DEFAULT_LIMITS;
  var state = {
    raw: "", chosenRecord: undefined, hosts: [S.defaultHost], goal: "BALANCED",
    natural: false, naturalOrganism: "",
    startChoice: "skip", startPaste: "", vectorAtg: false,
    nTagKey: "none", nTagCustom: "", cTagKey: "none", cTagCustom: "",
    presets: { biobrick: false, common: false }, ownSites: [],
    regions: [], regionsEdited: false,
    temperature: "37", temperatureEdited: false, hedge: false, seed: "",
    limits: { cai: "", mfe: "", gcMin: "", gcMax: "", stem: "" },
    advOpen: false, overrides: {}, seedUsed: null, running: false,
    lastReport: null, lastRequests: null, activeTab: 0, startParts: [],
  };
  var an = L.analyzeProtein("", {});
  var refs = {};

  function store(key, val) { try { sessionStorage.setItem(key, val); } catch (e) { /* private mode */ } }
  function load(key) { try { return sessionStorage.getItem(key); } catch (e) { return null; } }

  /* ---------- glossary tooltips ---------- */
  var GLOSS = {};
  S.glossary.terms.forEach(function (t) { GLOSS[t[0]] = t; });
  var openTip = null;
  function closeTip() { if (openTip) { openTip.btn.setAttribute("aria-expanded", "false"); openTip.bubble.hidden = true; openTip = null; } }
  document.addEventListener("keydown", function (e) { if (e.key === "Escape" && openTip) { var b = openTip.btn; closeTip(); b.focus(); } });
  document.addEventListener("click", function (e) { if (openTip && !openTip.wrap.contains(e.target)) closeTip(); });
  function tip(key) {
    var t = GLOSS[key];
    var bubbleId = id("tip");
    var bubble = h("span", { class: "bubble", role: "tooltip", id: bubbleId, hidden: true, text: t[2] });
    var btn = h("button", { type: "button", class: "tip", "aria-expanded": "false", "aria-describedby": bubbleId,
      "aria-label": F(S.glossary.help, { term: t[1] }), text: "?" });
    var wrap = h("span", { class: "tipwrap" }, btn, bubble);
    var pinned = false;
    function show() { if (openTip && openTip.btn !== btn) closeTip(); btn.setAttribute("aria-expanded", "true"); bubble.hidden = false; openTip = { btn: btn, bubble: bubble, wrap: wrap }; }
    function hide() { pinned = false; closeTip(); }
    // Click pins it open (or closes it); hover only previews; keyboard focus opens it.
    btn.addEventListener("click", function () { if (bubble.hidden) { pinned = true; show(); } else if (pinned) hide(); else pinned = true; });
    btn.addEventListener("focus", function () { if (btn.matches(":focus-visible")) show(); });
    wrap.addEventListener("mouseenter", show);
    wrap.addEventListener("mouseleave", function () { if (!pinned && !btn.matches(":focus-visible")) hide(); });
    btn.addEventListener("blur", function () { setTimeout(function () { if (openTip && openTip.btn === btn && !wrap.matches(":hover")) hide(); }, 100); });
    return wrap;
  }
  function labelWith(text, tips, forId) {
    var l = h("label", { class: "label", for: forId }, text);
    (tips || []).forEach(function (k) { l.appendChild(tip(k)); });
    return l;
  }

  /* ---------- live region ---------- */
  function announce(text) { refs.live.textContent = ""; setTimeout(function () { refs.live.textContent = text; }, 30); }

  /* ---------- clipboard ---------- */
  function copyText(text) {
    if (navigator.clipboard && navigator.clipboard.writeText) return navigator.clipboard.writeText(text);
    return new Promise(function (resolve, reject) {
      var ta = h("textarea", { "aria-hidden": "true", style: "position:fixed;opacity:0" }); ta.value = text;
      document.body.appendChild(ta); ta.select();
      try { document.execCommand("copy") ? resolve() : reject(); } catch (e) { reject(e); } finally { document.body.removeChild(ta); }
    });
  }

  function copyWithFeedback(text, msgEl, okText) {
    return copyText(text).then(function () { msgEl.textContent = okText; setTimeout(function () { msgEl.textContent = ""; }, 2500); })
      .catch(function () { msgEl.textContent = S.results.copyFailed; });
  }

  /* =====================================================================
     FORM
     ===================================================================== */
  function buildForm(seen) {
    var root = refs.form = h("div", { id: "formView" });
    if (!seen) root.appendChild(h("p", { class: "help", style: "font-size:1rem" }, S.app.intro));

    // ---- Step 1
    refs.ta = h("textarea", { id: "protein", spellcheck: "false", autocomplete: "off", autocapitalize: "characters",
      "aria-describedby": "protein-help protein-count protein-msgs" });
    refs.ta.addEventListener("input", function () { state.raw = refs.ta.value; state.chosenRecord = undefined; onProteinChange(); });
    refs.counter = h("span", { id: "protein-count", class: "muted small" });
    refs.msgs = h("div", { id: "protein-msgs" });
    refs.preview = h("div", { class: "preview", "aria-hidden": "true", hidden: true });
    var step1 = h("div", { class: "panel" },
      h("div", { class: "field", style: "margin-top:0" },
        labelWith(S.steps.protein.label, ["fasta"], "protein"),
        h("p", { class: "help", id: "protein-help" }, S.steps.protein.help),
        refs.ta,
        h("div", { class: "row between" },
          h("div", { class: "row", style: "margin:0" },
            h("button", { type: "button", onclick: function () { setProtein(S.app.exampleProtein); refs.ta.focus(); } }, S.app.tryExample),
            h("button", { type: "button", onclick: function () { setProtein(""); refs.ta.focus(); } }, S.app.clear)),
          refs.counter),
        refs.preview, refs.msgs));

    // ---- Step 2
    refs.hostCards = {};
    var commonGrid = h("div", { class: "cards four" }), moreGrid = h("div", { class: "cards" });
    HOST_IDS.forEach(function (hid) {
      var x = S.hosts[hid];
      var input = h("input", { type: "checkbox", value: hid, "aria-describedby": "d-" + hid });
      if (state.hosts.indexOf(hid) !== -1) input.checked = true;
      input.addEventListener("change", function () { onHostsChange(); });
      var card = h("label", { class: "card" }, input,
        h("span", { class: "body" },
          h("span", { class: "t" }, x.title),
          h("span", { class: "latin" }, x.latin),
          h("span", { class: "d", id: "d-" + hid }, x.desc)));
      refs.hostCards[hid] = { card: card, input: input, text: (x.title + " " + x.latin + " " + x.desc).toLowerCase() };
      (x.common ? commonGrid : moreGrid).appendChild(card);
    });
    refs.search = h("input", { type: "search", id: "host-search", "aria-describedby": "host-search-help" });
    refs.search.addEventListener("input", filterHosts);
    refs.noMatch = h("p", { class: "muted", hidden: true, role: "status" }, S.steps.host.noMatch);
    refs.moreDetails = h("details", { class: "more" },
      h("summary", null, S.steps.host.more),
      h("div", { class: "field" }, h("label", { for: "host-search", class: "label" }, S.steps.host.search),
        h("p", { class: "help", id: "host-search-help" }, S.steps.host.searchHelp), refs.search),
      refs.noMatch, moreGrid);
    refs.hostError = h("div", { class: "field-error", role: "alert", hidden: true });
    refs.tabsNote = h("p", { class: "help", hidden: true }, S.steps.host.tabsNote);
    var step2 = h("fieldset", { class: "panel" },
      h("legend", { class: "label" }, S.steps.host.label, tip("host")),
      h("p", { class: "help" }, S.steps.host.help), commonGrid, refs.moreDetails, refs.tabsNote, refs.hostError);

    // ---- Step 3
    var goalCards = h("div", { class: "cards", role: "radiogroup" });
    refs.goalInputs = {};
    ["PRODUCTION", "FOLDING", "BALANCED", "ALL"].forEach(function (g) {
      var x = S.goals[g];
      var input = h("input", { type: "radio", name: "goal", value: g, "aria-describedby": "gd-" + g });
      if (state.goal === g) input.checked = true;
      input.addEventListener("change", function () { state.goal = g; onGoalChange(); });
      refs.goalInputs[g] = input;
      goalCards.appendChild(h("label", { class: "card radio" }, input,
        h("span", { class: "body" },
          h("span", { class: "t" }, x.title, x.recommended ? h("span", { class: "rec" }, x.recommended) : null,
            h("span", { class: "tag", "aria-hidden": "true" }, x.tag)),
          h("span", { class: "d", id: "gd-" + g }, x.desc))));
    });
    refs.naturalBox = h("div", { class: "field" });
    var noR = h("input", { type: "radio", name: "natural", value: "no", checked: true });
    var yesR = h("input", { type: "radio", name: "natural", value: "yes" });
    refs.naturalInput = h("input", { type: "text", id: "natural-org", list: "org-list", autocomplete: "off", disabled: true });
    var dl = h("datalist", { id: "org-list" });
    HOST_IDS.forEach(function (hid) { dl.appendChild(h("option", { value: S.hosts[hid].title })); });
    refs.naturalError = h("div", { class: "field-error", hidden: true });
    refs.naturalOrgBox = h("div", { hidden: true, class: "field" },
      h("label", { class: "label", for: "natural-org" }, S.steps.natural.organism),
      h("p", { class: "help" }, S.steps.natural.organismHelp), refs.naturalInput, dl, refs.naturalError);
    function onNat() {
      state.natural = yesR.checked; refs.naturalOrgBox.hidden = !state.natural; refs.naturalInput.disabled = !state.natural;
      if (!state.natural) { state.naturalOrganism = ""; refs.naturalInput.value = ""; refs.naturalError.hidden = true; }
    }
    noR.addEventListener("change", onNat); yesR.addEventListener("change", onNat);
    refs.naturalInput.addEventListener("input", function () {
      var v = refs.naturalInput.value.trim().toLowerCase();
      var hit = HOST_IDS.filter(function (hid) { return S.hosts[hid].title.toLowerCase() === v || S.hosts[hid].latin.toLowerCase() === v; })[0];
      state.naturalOrganism = hit || ""; refs.naturalError.hidden = true;
    });
    refs.naturalBox.appendChild(h("fieldset", null,
      h("legend", { class: "label" }, S.steps.natural.label, tip("denovo")),
      h("div", { class: "row" },
        h("label", { class: "check" }, noR, S.steps.natural.no),
        h("label", { class: "check" }, yesR, S.steps.natural.yes)),
      refs.naturalOrgBox));
    var step3 = h("fieldset", { class: "panel" }, h("legend", { class: "label" }, S.steps.goal.label), goalCards, refs.naturalBox);

    // ---- Action bar
    refs.optimize = h("button", { type: "button", class: "primary", id: "optimize", "aria-disabled": "true", title: S.app.optimizeDisabledTip, "aria-describedby": "optimize-note" }, S.app.optimize);
    refs.optimize.addEventListener("click", onOptimizeClick);
    refs.advBadge = h("button", { type: "button", class: "link", "aria-expanded": "false", "aria-controls": "adv-panel" }, S.app.advancedLink);
    refs.advBadge.addEventListener("click", function () { setAdvanced(!state.advOpen); });
    refs.formError = h("div", { class: "msg error", role: "alert", hidden: true });
    refs.actions = h("div", { class: "actions-bar" }, refs.optimize,
      h("p", { class: "help", id: "optimize-note" }, S.app.belowButton, " ", h("a", { href: "/glossary" }, S.app.glossaryLink)),
      refs.formError);
    refs.loading = h("div", { class: "panel loading", hidden: true, role: "status" });

    var adv = h("div", { class: "adv-toggle" }, refs.advBadge);
    refs.advPanel = h("div", { class: "panel adv-panel", id: "adv-panel", hidden: true });
    buildAdvanced(refs.advPanel);

    add(root, [step1, step2, step3, refs.actions, refs.loading, adv, refs.advPanel]);
    return root;
  }

  function setProtein(text) { refs.ta.value = text; state.raw = text; state.chosenRecord = undefined; onProteinChange(); }

  /* ---------- step 1 updates ---------- */
  function msgEl(kind, ico, text, actions) {
    return h("div", { class: "msg " + kind }, h("span", { class: "ico", "aria-hidden": "true" }, ico), h("span", null, text),
      actions && actions.length ? h("span", { class: "actions" }, actions) : null);
  }
  function clearFormError() { if (refs.formError) refs.formError.hidden = true; }
  function onProteinChange() {
    clearFormError();
    an = L.analyzeProtein(state.raw, { maxLength: api.max_protein_length, chosenRecord: state.chosenRecord });
    clear(refs.msgs);
    var n = an.sequence.replace(/\*/g, "").length;
    refs.counter.textContent = an.empty || an.needsRecordChoice ? "" : (n === 1 ? S.steps.protein.counterOne : F(S.steps.protein.counter, { n: n.toLocaleString("en-US") }));
    refs.ta.setAttribute("aria-invalid", an.errors.length && !an.empty ? "true" : "false");

    an.errors.forEach(function (e) {
      var actions = [];
      if (e.id === "invalidLetters") actions.push(h("button", { type: "button", onclick: function () { setProtein(L.removeInvalid(an.sequence)); } }, S.errors.removeThem));
      if (e.id === "dna") actions.push(h("button", { type: "button", onclick: function () { setProtein(L.translateFrame1(an.dnaSequence)); } }, S.errors.translate));
      if (e.id === "multi") {
        actions.push(h("button", { type: "button", onclick: function () { state.chosenRecord = 0; onProteinChange(); } }, S.errors.useFirst));
        var sel = h("select", { "aria-label": S.errors.chooseLabel });
        sel.appendChild(h("option", { value: "" }, S.errors.chooseOne));
        an.recordNames.forEach(function (nm, i) { sel.appendChild(h("option", { value: String(i) }, nm)); });
        sel.addEventListener("change", function () { if (sel.value !== "") { state.chosenRecord = Number(sel.value); onProteinChange(); } });
        actions.push(sel);
      }
      refs.msgs.appendChild(msgEl("error", "✖", e.text, actions));
    });
    an.warnings.forEach(function (w) { refs.msgs.appendChild(msgEl("warn", "⚠", w.text)); });
    an.infos.forEach(function (w) { refs.msgs.appendChild(msgEl("info", "ℹ", w.text)); });
    an.notes.forEach(function (w) { refs.msgs.appendChild(msgEl("note", "✓", w.text)); });

    if (an.invalid.length) {
      var bad = {}; an.invalid.forEach(function (x) { bad[x.index] = true; });
      clear(refs.preview); var run = "";
      for (var i = 0; i < an.sequence.length; i++) {
        if (bad[i]) { if (run) { refs.preview.appendChild(document.createTextNode(run)); run = ""; } refs.preview.appendChild(h("mark", null, an.sequence[i])); }
        else run += an.sequence[i];
      }
      if (run) refs.preview.appendChild(document.createTextNode(run));
      refs.preview.hidden = false;
    } else refs.preview.hidden = true;

    if (an.errors.length && !an.empty) announce(an.errors.length === 1 ? S.errors.problemsOne : F(S.errors.problems, { n: an.errors.length }));

    if (!state.regionsEdited) { state.regions = L.detectRegions(an.sequence); renderRegions(); }
    updateOptimizeState(); updateBadge();
  }
  function updateOptimizeState() {
    var ok = !an.blocked && !an.empty;
    refs.optimize.setAttribute("aria-disabled", ok ? "false" : "true");
    if (ok) refs.optimize.removeAttribute("title"); else refs.optimize.setAttribute("title", S.app.optimizeDisabledTip);
  }

  /* ---------- step 2 / 3 updates ---------- */
  function onHostsChange() {
    clearFormError();
    state.hosts = HOST_IDS.filter(function (hid) { return refs.hostCards[hid].input.checked; });
    refs.hostError.hidden = state.hosts.length > 0; if (!state.hosts.length) refs.hostError.textContent = S.errors.noHost;
    refs.tabsNote.hidden = state.hosts.length <= 3;
    syncTemperature(); updateBadge(); refreshStartParts();
  }
  function filterHosts() {
    var q = refs.search.value.trim().toLowerCase(), any = false;
    HOST_IDS.forEach(function (hid) { var m = !q || refs.hostCards[hid].text.indexOf(q) !== -1; refs.hostCards[hid].card.hidden = !m; any = any || m; });
    refs.noMatch.hidden = any;
    if (q) refs.moreDetails.open = true;
  }
  function onGoalChange() {
    refs.naturalBox.hidden = state.goal === "PRODUCTION";
    if (state.goal === "PRODUCTION") { state.natural = false; state.naturalOrganism = ""; }
    renderRegions();
  }

  /* =====================================================================
     ADVANCED
     ===================================================================== */
  function buildAdvanced(panel) {
    var A = S.advanced;

    // ---- Your vector
    refs.startSelect = h("select", { id: "start-choice" });
    refs.startSelect.appendChild(h("option", { value: "skip" }, A.start.skip));
    refs.startSelect.appendChild(h("option", { value: "paste" }, A.start.paste));
    refs.startGroup = h("optgroup", { label: A.start.partsGroup });
    refs.startSelect.appendChild(refs.startGroup);
    refs.startPaste = h("input", { type: "text", id: "start-paste", autocomplete: "off", spellcheck: "false", "aria-describedby": "start-paste-err" });
    refs.startPasteErr = h("div", { class: "field-error", id: "start-paste-err", hidden: true });
    refs.startPasteNote = h("div", { class: "msg note", hidden: true });
    refs.startPasteBox = h("div", { class: "field", hidden: true }, h("label", { class: "label", for: "start-paste" }, A.start.pasteLabel), refs.startPaste, refs.startPasteErr, refs.startPasteNote);
    refs.startSelect.addEventListener("change", function () {
      state.startChoice = refs.startSelect.value; refs.startPasteBox.hidden = state.startChoice !== "paste"; validateStart(); updateBadge();
    });
    refs.startPaste.addEventListener("input", function () { state.startPaste = refs.startPaste.value; updateBadge(); });
    refs.startPaste.addEventListener("blur", function () { validateStart(); updateBadge(); });
    refs.vectorAtg = h("input", { type: "checkbox", id: "vector-atg" });
    refs.vectorAtg.addEventListener("change", function () { state.vectorAtg = refs.vectorAtg.checked; updateBadge(); });
    var gVector = group(A.groups.vector,
      h("div", { class: "field", style: "margin-top:0" }, labelWith(A.start.label, ["rbs"], "start-choice"), h("p", { class: "help" }, A.start.help), refs.startSelect),
      refs.startPasteBox,
      h("label", { class: "check" }, refs.vectorAtg, h("span", null, A.vectorAtg.label, tip("atg"), h("span", { class: "help", style: "display:block;margin:0" }, A.vectorAtg.help))));

    // ---- Tags
    function tagField(which, spec) {
      var selId = id("tag");
      var sel = h("select", { id: selId });
      Object.keys(A.tags).forEach(function (k) { sel.appendChild(h("option", { value: k }, A.tags[k].name + (A.tags[k].aa ? " — " + A.tags[k].aa : ""))); });
      var custom = h("input", { type: "text", id: selId + "c", autocomplete: "off", spellcheck: "false", "aria-describedby": selId + "e" });
      var err = h("div", { class: "field-error", id: selId + "e", hidden: true });
      var box = h("div", { class: "field", hidden: true }, h("label", { class: "label", for: selId + "c" }, A.tagCustomLabel), custom, err);
      refs[which] = { sel: sel, custom: custom, err: err, box: box };
      sel.addEventListener("change", function () { state[which + "Key"] = sel.value; box.hidden = sel.value !== "custom"; validateTags(); updateBadge(); });
      custom.addEventListener("input", function () { state[which + "Custom"] = custom.value; updateBadge(); });
      custom.addEventListener("blur", function () { validateTags(); updateBadge(); });
      return h("div", { class: "field" }, labelWith(spec.label, which === "nTag" ? ["tag", "tev"] : ["tag"], selId), h("p", { class: "help" }, spec.help), sel, box);
    }
    var gTags = group(A.groups.tags, tagField("nTag", A.nTag), tagField("cTag", A.cTag));

    // ---- Cut sites
    var presetBox = function (key, text) {
      var c = h("input", { type: "checkbox" });
      c.addEventListener("change", function () { state.presets[key] = c.checked; updateBadge(); });
      refs["preset_" + key] = c;
      return h("label", { class: "check" }, c, text);
    };
    refs.chips = h("div", { class: "chips" });
    refs.siteName = h("input", { type: "text", id: "site-name", list: "enzyme-list", autocomplete: "off" });
    var enzDl = h("datalist", { id: "enzyme-list" });
    Object.keys(A.enzymes).forEach(function (n) { enzDl.appendChild(h("option", { value: n }, A.enzymes[n])); });
    refs.siteSeq = h("input", { type: "text", id: "site-seq", autocomplete: "off", spellcheck: "false", "aria-describedby": "site-err" });
    refs.siteErr = h("div", { class: "field-error", id: "site-err", hidden: true });
    refs.siteName.addEventListener("change", function () { var k = refs.siteName.value.trim(); if (A.enzymes[k] && !refs.siteSeq.value) refs.siteSeq.value = A.enzymes[k]; });
    var addBtn = h("button", { type: "button", onclick: addSite }, A.sites.add);
    var gSites = group(A.groups.sites,
      h("div", { class: "field", style: "margin-top:0" }, h("span", { class: "label" }, A.sites.label, tip("cutsite"), tip("goldengate"), tip("biobrick")), h("p", { class: "help" }, A.sites.help)),
      h("p", { class: "msg note" }, h("span", { class: "ico", "aria-hidden": "true" }, "🔒"), A.sites.locked),
      presetBox("biobrick", A.sites.biobrick), presetBox("common", A.sites.common),
      h("div", { class: "field" }, h("span", { class: "label" }, A.sites.addOwn),
        h("div", { class: "row" },
          h("div", { style: "flex:1;min-width:140px" }, h("label", { for: "site-name", class: "help" }, A.sites.addName), refs.siteName, enzDl),
          h("div", { style: "flex:1;min-width:140px" }, h("label", { for: "site-seq", class: "help" }, A.sites.addSite), refs.siteSeq),
          h("div", { style: "align-self:flex-end" }, addBtn)),
        refs.siteErr, refs.chips));

    // ---- Regions
    refs.regionSeq = h("div", { class: "preview", hidden: true });
    refs.regionBody = h("tbody"); refs.regionNone = h("p", { class: "muted" }, A.regions.none);
    refs.regionRelevant = h("p", { class: "help", hidden: true }, A.regions.relevant);
    var regionTable = h("table", { class: "regions stack" },
      h("thead", null, h("tr", null, h("th", null, A.regions.start), h("th", null, A.regions.end), h("th", null, A.regions.type), h("th", null, ""))), refs.regionBody);
    var gRegions = group(A.groups.regions,
      h("div", { class: "field", style: "margin-top:0" }, h("span", { class: "label" }, A.regions.label, tip("linker"), tip("folding")), h("p", { class: "help" }, A.regions.help), refs.regionRelevant),
      refs.regionSeq, regionTable, refs.regionNone,
      h("div", { class: "row" },
        h("button", { type: "button", onclick: function () { state.regions.push({ start: 1, end: 1, type: "linker" }); state.regionsEdited = true; renderRegions(); updateBadge(); } }, A.regions.add),
        h("button", { type: "button", onclick: function () { state.regionsEdited = false; state.regions = L.detectRegions(an.sequence); renderRegions(); updateBadge(); } }, A.regions.reset)));

    // ---- Growth
    refs.temp = h("input", { type: "number", id: "temp", step: "0.5", inputmode: "decimal", "aria-describedby": "temp-help temp-err" });
    refs.tempErr = h("div", { class: "field-error", id: "temp-err", hidden: true });
    refs.tempNote = h("p", { class: "help", id: "temp-help" });
    refs.quick = h("div", { class: "quick", role: "group", "aria-label": A.temperature.quick });
    [37, 30, 25, 18].forEach(function (t) {
      var b = h("button", { type: "button", "aria-pressed": "false", "data-t": String(t) }, t + " °C");
      b.addEventListener("click", function () { refs.temp.value = String(t); onTempInput(); });
      refs.quick.appendChild(b);
    });
    refs.temp.addEventListener("input", onTempInput);
    var gGrowth = group(A.groups.growth,
      h("div", { class: "field", style: "margin-top:0" }, h("label", { class: "label", for: "temp" }, A.temperature.label), h("p", { class: "help" }, A.temperature.help),
        refs.temp, refs.tempErr, refs.tempNote, refs.quick));

    // ---- Results
    refs.hedge = h("input", { type: "checkbox", id: "hedge" });
    refs.hedge.addEventListener("change", function () { state.hedge = refs.hedge.checked; updateBadge(); });
    refs.seed = h("input", { type: "text", id: "seed", inputmode: "numeric", autocomplete: "off", "aria-describedby": "seed-err" });
    refs.seedErr = h("div", { class: "field-error", id: "seed-err", hidden: true });
    refs.seed.addEventListener("input", function () { state.seed = refs.seed.value; updateBadge(); });
    refs.seed.addEventListener("blur", function () { validateSeed(); updateBadge(); });
    var gResults = group(A.groups.results,
      h("label", { class: "check" }, refs.hedge, h("span", null, A.backup.label, tip("backup"), h("span", { class: "help", style: "display:block;margin:0" }, A.backup.help))),
      h("div", { class: "field" }, labelWith(A.seed.label, ["seed"], "seed"), h("p", { class: "help" }, A.seed.help), refs.seed, refs.seedErr));

    // ---- Expert limits
    refs.lim = {};
    function limField(key, label, terms, step) {
      var i = h("input", { type: "number", id: "lim-" + key, step: step, inputmode: "decimal" });
      refs.lim[key] = i;
      i.addEventListener("input", function () { state.limits[key] = i.value; refs.limWarn.hidden = !limitsChanged(); updateBadge(); });
      return h("div", { class: "field", style: "flex:1;min-width:200px" }, labelWith(label, terms, "lim-" + key), i);
    }
    refs.limWarn = h("p", { class: "msg warn", role: "status", hidden: true }, h("span", { class: "ico", "aria-hidden": "true" }, "⚠"), A.limits.warning);
    var gLimits = group(A.groups.limits,
      h("p", { class: "help", style: "margin-top:0" }, A.limits.warning),
      h("div", { class: "row" }, limField("cai", A.limits.cai, ["cai"], "0.01"), limField("mfe", A.limits.mfe, ["hairpin", "fold"], "1")),
      h("div", { class: "row" }, limField("gcMin", A.limits.gcMin, ["gc"], "1"), limField("gcMax", A.limits.gcMax, ["gc"], "1"), limField("stem", A.limits.stem, ["hairpin"], "1")),
      refs.limWarn);

    // ---- Developer options
    refs.regionsJson = h("textarea", { id: "regions-json", style: "min-height:90px", spellcheck: "false" });
    refs.jsonMsg = h("div", { class: "msg", hidden: true });
    refs.copyMsg = h("span", { class: "muted small", role: "status" });
    var dev = h("details", { class: "adv-group" }, h("summary", { class: "label" }, A.developer),
      h("div", { class: "row" },
        h("button", { type: "button", onclick: exportRegions }, A.developerOptions.export),
        h("button", { type: "button", onclick: importRegions }, A.developerOptions.importApply),
        h("button", { type: "button", onclick: copyApiRequest }, A.developerOptions.copyRequest), refs.copyMsg),
      h("div", { class: "field" }, h("label", { class: "label", for: "regions-json" }, A.developerOptions.importLabel), refs.regionsJson, refs.jsonMsg));

    refs.reset = h("button", { type: "button", class: "link", onclick: resetAdvanced }, A.reset);
    add(panel, [h("div", { class: "row between", style: "margin-top:0" }, h("h2", null, S.app.advancedLink), refs.reset),
      gVector, gTags, gSites, gRegions, gGrowth, gResults, gLimits, dev]);
    syncTemperature();
  }
  function group(title, /* children */) {
    var g = h("div", { class: "adv-group" }, h("h3", null, title));
    for (var i = 1; i < arguments.length; i++) add(g, arguments[i]);
    return g;
  }

  /* ---------- advanced logic ---------- */
  function setAdvanced(open) {
    state.advOpen = open; refs.advPanel.hidden = !open; refs.advBadge.setAttribute("aria-expanded", open ? "true" : "false");
    store("pd_adv", open ? "1" : "0");
  }
  function limitsChanged() {
    var d = api.defaults;
    var map = { cai: d.cai_floor, mfe: d.worst_window_energy, gcMin: d.gc_min * 100, gcMax: d.gc_max * 100, stem: d.longest_allowed_stem };
    return Object.keys(map).some(function (k) { var v = state.limits[k]; return v !== "" && v !== null && Number(v) !== map[k]; });
  }
  function tagValue(which) {
    var key = state[which + "Key"];
    return key === "custom" ? state[which + "Custom"] : (S.advanced.tags[key] || {}).aa || "";
  }
  function changedCount() {
    var n = 0;
    if (state.startChoice !== "skip") n++;
    if (state.vectorAtg) n++;
    if (state.nTagKey !== "none") n++;
    if (state.cTagKey !== "none") n++;
    if (state.presets.biobrick || state.presets.common || state.ownSites.length) n++;
    if (state.regionsEdited) n++;
    if (state.temperatureEdited) n++;
    if (state.hedge) n++;
    if (String(state.seed).trim() !== "") n++;
    if (limitsChanged()) n++;
    return n;
  }
  function updateBadge() {
    var n = changedCount();
    refs.advBadge.textContent = n ? F(S.app.advancedChanged, { n: n }) : S.app.advancedLink;
  }
  function suggestedTemps() {
    return state.hosts.map(function (hid) { return L.hostTemperature(hid); });
  }
  function showLimitDefaults() {
    var d = api.defaults, map = { cai: d.cai_floor.toFixed(2), mfe: String(d.worst_window_energy), gcMin: String(Math.round(d.gc_min * 100)),
      gcMax: String(Math.round(d.gc_max * 100)), stem: String(d.longest_allowed_stem) };
    Object.keys(map).forEach(function (k) { if (refs.lim && refs.lim[k]) refs.lim[k].setAttribute("placeholder", map[k]); });
  }
  function refreshStartParts() {
    var host = state.hosts[0];
    if (!host || !refs.startGroup) return;
    fetch("/api/rbs-options?host=" + encodeURIComponent(host)).then(function (r) { return r.json(); }).then(function (d) {
      if (state.hosts[0] !== host) return; // the user picked another host meanwhile
      state.startParts = d.parts || [];
      clear(refs.startGroup);
      state.startParts.forEach(function (p) { refs.startGroup.appendChild(h("option", { value: "part:" + p.name }, p.name + " \u2014 " + p.dna)); });
      if (!state.startParts.length) refs.startGroup.appendChild(h("option", { value: "", disabled: true }, S.advanced.start.partsNone));
      if (state.startChoice.indexOf("part:") === 0 && !state.startParts.some(function (p) { return "part:" + p.name === state.startChoice; })) {
        state.startChoice = "skip"; refs.startSelect.value = "skip"; updateBadge();
      } else refs.startSelect.value = state.startChoice;
    }).catch(function () { /* keep the current list */ });
  }
  function syncTemperature() {
    if (!refs.temp) return;
    var temps = suggestedTemps(), mixed = temps.some(function (t) { return t !== temps[0]; });
    if (!state.temperatureEdited) { state.temperature = String(temps.length ? temps[0] : 37); refs.temp.value = state.temperature; }
    refs.tempNote.textContent = state.temperatureEdited ? "" : (mixed ? S.advanced.temperature.mixed : S.advanced.temperature.hostNote);
    markQuick();
  }
  function markQuick() {
    Array.prototype.forEach.call(refs.quick.children, function (b) { b.setAttribute("aria-pressed", String(Number(refs.temp.value) === Number(b.getAttribute("data-t")))); });
  }
  function onTempInput() {
    state.temperature = refs.temp.value; state.temperatureEdited = true; markQuick();
    var v = L.validateTemperature(state.temperature);
    refs.tempErr.hidden = v.ok; refs.tempErr.textContent = v.ok ? "" : v.error; refs.temp.setAttribute("aria-invalid", v.ok ? "false" : "true");
    refs.tempNote.textContent = ""; updateBadge();
  }
  function validateSeed() {
    var v = L.validateSeed(state.seed);
    refs.seedErr.hidden = v.ok; refs.seedErr.textContent = v.ok ? "" : v.error; refs.seed.setAttribute("aria-invalid", v.ok ? "false" : "true");
    return v.ok;
  }
  function validateStart() {
    if (state.startChoice !== "paste") { refs.startPasteErr.hidden = true; refs.startPasteNote.hidden = true; return true; }
    var v = L.validateStartDna(state.startPaste);
    refs.startPasteErr.hidden = v.ok; refs.startPasteErr.textContent = v.ok ? "" : v.error;
    refs.startPaste.setAttribute("aria-invalid", v.ok ? "false" : "true");
    refs.startPasteNote.hidden = !(v.ok && v.note); refs.startPasteNote.textContent = v.note || "";
    return v.ok;
  }
  function validateTags() {
    var ok = true;
    ["nTag", "cTag"].forEach(function (w) {
      var r = refs[w], key = state[w + "Key"];
      if (key !== "custom") { r.err.hidden = true; return; }
      var v = L.validateTag(state[w + "Custom"]);
      r.err.hidden = v.ok && !v.warning; r.err.textContent = v.ok ? (v.warning || "") : v.error;
      r.custom.setAttribute("aria-invalid", v.ok ? "false" : "true"); if (!v.ok) ok = false;
    });
    return ok;
  }
  function resetAdvanced() {
    state.startChoice = "skip"; state.startPaste = ""; state.vectorAtg = false; state.nTagKey = "none"; state.cTagKey = "none";
    state.nTagCustom = ""; state.cTagCustom = ""; state.presets = { biobrick: false, common: false }; state.ownSites = [];
    state.regionsEdited = false; state.regions = L.detectRegions(an.sequence); state.temperatureEdited = false;
    state.hedge = false; state.seed = ""; state.limits = { cai: "", mfe: "", gcMin: "", gcMax: "", stem: "" }; state.overrides = {};
    refs.startSelect.value = "skip"; refs.startPaste.value = ""; refs.startPasteBox.hidden = true; refs.vectorAtg.checked = false;
    ["nTag", "cTag"].forEach(function (w) { refs[w].sel.value = "none"; refs[w].custom.value = ""; refs[w].box.hidden = true; refs[w].err.hidden = true; });
    refs.preset_biobrick.checked = false; refs.preset_common.checked = false; renderChips();
    refs.hedge.checked = false; refs.seed.value = ""; refs.seedErr.hidden = true; refs.tempErr.hidden = true;
    Object.keys(refs.lim).forEach(function (k) { refs.lim[k].value = ""; }); refs.limWarn.hidden = true;
    refs.startPasteErr.hidden = true; refs.startPasteNote.hidden = true;
    syncTemperature(); renderRegions(); updateBadge();
  }

  /* ---------- cut sites ---------- */
  function addSite() {
    var name = refs.siteName.value.trim(), site = refs.siteSeq.value.trim();
    if (!site && S.advanced.enzymes[name]) site = S.advanced.enzymes[name];
    var v = L.validateSite(site);
    if (!v.ok) { refs.siteErr.hidden = false; refs.siteErr.textContent = v.error; refs.siteSeq.setAttribute("aria-invalid", "true"); return; }
    if (!name) { refs.siteErr.hidden = false; refs.siteErr.textContent = S.errors.siteName; return; }
    refs.siteErr.hidden = true; refs.siteSeq.removeAttribute("aria-invalid");
    state.ownSites.push({ name: name, site: v.site });
    refs.siteName.value = ""; refs.siteSeq.value = ""; renderChips(); updateBadge();
  }
  function renderChips() {
    clear(refs.chips);
    state.ownSites.forEach(function (s, i) {
      refs.chips.appendChild(h("span", { class: "chip" }, s.name + " " + s.site,
        h("button", { type: "button", "aria-label": F(S.advanced.sites.remove, { name: s.name }), onclick: function () { state.ownSites.splice(i, 1); renderChips(); updateBadge(); refs.siteName.focus(); } }, "×")));
    });
  }
  function allSites() {
    var out = [], seen = {};
    function push(n, s) { if (!seen[n + s]) { seen[n + s] = 1; out.push({ name: n, site: s }); } }
    ["biobrick", "common"].forEach(function (k) { if (state.presets[k]) S.advanced.presets[k].forEach(function (n) { push(n, S.advanced.enzymes[n]); }); });
    state.ownSites.forEach(function (s) { push(s.name, s.site); });
    return out;
  }

  /* ---------- regions ---------- */
  function renderRegions() {
    clear(refs.regionBody);
    var errs = L.validateRegions(state.regions, an.sequence.length || 0);
    state.regions.forEach(function (r, i) {
      function num(key) {
        var el = h("input", { type: "number", min: "1", value: String(r[key]), "aria-label": S.advanced.regions[key] + " " + (i + 1) });
        el.addEventListener("input", function () { r[key] = Number(el.value); state.regionsEdited = true; updateRegionErrors(); updateBadge(); });
        return el;
      }
      var sel = h("select", { "aria-label": S.advanced.regions.type + " " + (i + 1) });
      ["linker", "loop", "domain_boundary"].forEach(function (k) { var o = h("option", { value: k }, S.advanced.regions[k]); if (r.type === k) o.selected = true; sel.appendChild(o); });
      sel.addEventListener("change", function () { r.type = sel.value; state.regionsEdited = true; updateBadge(); });
      var err = h("div", { class: "field-error", role: "alert", hidden: !errs[i] }, errs[i] || "");
      var tr = h("tr", { "data-row": String(i) },
        h("td", { "data-label": S.advanced.regions.start }, num("start")),
        h("td", { "data-label": S.advanced.regions.end }, num("end")),
        h("td", { "data-label": S.advanced.regions.type }, sel),
        h("td", null, h("button", { type: "button", "aria-label": F(S.errors.removeRegion, { n: i + 1 }), onclick: function () { state.regions.splice(i, 1); state.regionsEdited = true; renderRegions(); updateBadge(); } }, S.advanced.regions.remove)));
      refs.regionBody.appendChild(tr);
      refs.regionBody.appendChild(h("tr", { hidden: !errs[i], "data-err": String(i) }, h("td", { colspan: "4" }, err)));
    });
    refs.regionNone.hidden = state.regions.length > 0;
    refs.regionBody.parentNode.hidden = state.regions.length === 0;
    refs.regionRelevant.hidden = state.goal !== "PRODUCTION";
    // sequence with linkers highlighted
    var seq = an.sequence, auto = L.detectRegions(seq);
    if (seq && auto.length) {
      clear(refs.regionSeq); var pos = 0;
      auto.forEach(function (r) { refs.regionSeq.appendChild(document.createTextNode(seq.slice(pos, r.start - 1))); refs.regionSeq.appendChild(h("mark", { style: "background:var(--accent)" }, seq.slice(r.start - 1, r.end))); pos = r.end; });
      refs.regionSeq.appendChild(document.createTextNode(seq.slice(pos))); refs.regionSeq.setAttribute("aria-label", S.advanced.regions.highlighted); refs.regionSeq.hidden = false;
    } else refs.regionSeq.hidden = true;
  }
  function updateRegionErrors() {
    var errs = L.validateRegions(state.regions, an.sequence.length || 0);
    errs.forEach(function (e, i) {
      var row = refs.regionBody.querySelector('tr[data-err="' + i + '"]');
      if (!row) return; row.hidden = !e; row.querySelector(".field-error").textContent = e || ""; row.querySelector(".field-error").hidden = !e;
    });
  }
  function exportRegions() {
    refs.regionsJson.value = JSON.stringify(state.regions.map(function (r) { return { start: r.start, end: r.end, kind: r.type }; }), null, 2);
  }
  function importRegions() {
    var box = refs.jsonMsg; box.hidden = true;
    try {
      var data = JSON.parse(refs.regionsJson.value);
      if (!Array.isArray(data) || data.some(function (d) { return typeof d.start !== "number" || typeof d.end !== "number" || ["linker", "loop", "domain_boundary"].indexOf(d.kind) === -1; })) throw new Error("shape");
      state.regions = data.map(function (d) { return { start: d.start, end: d.end, type: d.kind }; }); state.regionsEdited = true; renderRegions(); updateBadge();
    } catch (e) { box.className = "msg error"; box.textContent = S.advanced.developerOptions.importBad; box.hidden = false; }
  }

  /* =====================================================================
     RUN
     ===================================================================== */
  function collectState() {
    var startDna = "";
    if (state.startChoice === "paste") startDna = L.validateStartDna(state.startPaste).dna || "";
    else if (state.startChoice.indexOf("part:") === 0) {
      var part = state.startParts.filter(function (p) { return "part:" + p.name === state.startChoice; })[0]; startDna = part ? part.dna : "";
    }
    var lim = {
      cai: state.limits.cai, mfe: state.limits.mfe, stem: state.limits.stem,
      gcMin: state.limits.gcMin === "" ? "" : Number(state.limits.gcMin) / 100,
      gcMax: state.limits.gcMax === "" ? "" : Number(state.limits.gcMax) / 100,
    };
    return {
      protein: an.sequence.replace(/\*/g, ""), hosts: state.hosts.slice(), goal: state.goal,
      naturalOrganism: state.natural ? state.naturalOrganism : "", startDna: startDna, vectorAtg: state.vectorAtg,
      nTag: L.validateTag(tagValue("nTag")).sequence, cTag: L.validateTag(tagValue("cTag")).sequence,
      sites: allSites(), regions: state.regions.filter(function (r) { return r.end >= r.start; }), hedge: state.hedge,
      seed: (function () { var v = L.validateSeed(state.seed); return v.ok && v.value !== null ? v.value : state.seedUsed; })(), temperature: state.temperature, temperatureEdited: state.temperatureEdited,
      limits: lim, overrides: state.overrides,
    };
  }
  function validateAll() {
    var problems = [], first = null;
    function bad(el) { problems.push(el); if (!first) first = el; }
    if (an.blocked || an.empty) bad(refs.ta);
    if (!state.hosts.length) { refs.hostError.hidden = false; refs.hostError.textContent = S.errors.noHost; bad(refs.hostCards[HOST_IDS[0]].input); }
    if (state.natural && !state.naturalOrganism) { refs.naturalError.hidden = false; refs.naturalError.textContent = S.errors.organismRequired; bad(refs.naturalInput); }
    if (!validateSeed()) bad(refs.seed);
    if (!validateStart()) bad(refs.startPaste);
    if (!validateTags()) bad(refs.nTag.custom);
    if (state.temperatureEdited && !L.validateTemperature(state.temperature).ok) bad(refs.temp);
    var rerr = L.validateRegions(state.regions, an.sequence.length).some(Boolean); if (rerr) { updateRegionErrors(); bad(refs.regionBody); }
    return { ok: problems.length === 0, first: first, count: problems.length };
  }
  function onOptimizeClick() {
    if (state.running) return;
    if (refs.optimize.getAttribute("aria-disabled") === "true") {
      clear(refs.formError); refs.formError.hidden = false;
      refs.formError.appendChild(h("span", { class: "ico", "aria-hidden": "true" }, "✖"));
      refs.formError.appendChild(document.createTextNode(an.empty ? S.errors.empty : S.errors.fixFields));
      refs.ta.focus();
      return;
    }
    var v = validateAll();
    if (!v.ok) {
      clear(refs.formError); refs.formError.hidden = false;
      refs.formError.appendChild(h("span", { class: "ico", "aria-hidden": "true" }, "✖"));
      refs.formError.appendChild(document.createTextNode(S.errors.fixFields));
      announce(v.count === 1 ? S.errors.problemsOne : F(S.errors.problems, { n: v.count }));
      if (v.first) { if (v.first === refs.hostCards[HOST_IDS[0]].input) v.first.closest("fieldset").scrollIntoView({ block: "center" }); if (v.first.focus) v.first.focus(); }
      return;
    }
    refs.formError.hidden = true;
    state.overrides = {};
    run(true);
  }

  var controller = null, timer = null;
  function run(fresh) {
    if (state.running) return;
    var seedCheck = L.validateSeed(state.seed);
    if (fresh || state.seedUsed === null) state.seedUsed = seedCheck.value !== null ? seedCheck.value : Math.floor(Math.random() * 1e6) + 1;
    var requests = L.buildRequests(collectState(), api.defaults);
    state.lastRequests = requests; state.running = true;
    showLoading(true);
    controller = new AbortController();
    Promise.all(requests.map(function (body) {
      return fetch("/api/optimize", { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(body), signal: controller.signal })
        .then(function (res) { return res.json().then(function (data) { if (!res.ok) { var e = new Error(data.error || "server"); e.userMessage = data.error; e.status = res.status; throw e; } return data; }); });
    })).then(function (reports) {
      state.lastReport = L.mergeReports(reports); state.activeTab = 0; showLoading(false);
      try { localStorage.setItem("pd_seen", "1"); } catch (e) { /* ignore */ }
      renderResults(); showView("results");
    }).catch(function (err) {
      showLoading(false);
      if (err && err.name === "AbortError") { announce(S.loading.cancelled); refs.optimize.focus(); return; }
      showRunError(err);
    });
  }
  function showLoading(on) {
    state.running = on; refs.optimize.hidden = on; refs.loading.hidden = !on;
    clearInterval(timer);
    if (!on) return;
    clear(refs.loading);
    var elapsed = h("p", { class: "muted" }), slow = h("p", { class: "msg warn", hidden: true }, S.loading.slow);
    var list = h("ol", { "aria-label": S.loading.title }); S.loading.steps.forEach(function (t) { list.appendChild(h("li", null, t)); });
    add(refs.loading, [h("h2", null, h("span", { class: "spinner", "aria-hidden": "true" }), S.loading.title), list, elapsed, slow,
      h("button", { type: "button", onclick: function () { if (controller) controller.abort(); } }, S.loading.cancel)]);
    var t0 = Date.now();
    function tick() { var s = Math.floor((Date.now() - t0) / 1000); elapsed.textContent = F(S.loading.elapsed, { s: s }); if (s >= 60) slow.hidden = false; }
    tick(); timer = setInterval(tick, 1000);
    refs.loading.setAttribute("tabindex", "-1"); refs.loading.focus();
    announce(S.loading.title);
  }
  function showRunError(err) {
    clear(refs.formError); refs.formError.hidden = false;
    var offline = typeof navigator !== "undefined" && navigator.onLine === false;
    var text = offline ? S.errors.offline : (err && err.status === 400 && err.userMessage ? err.userMessage : S.errors.server);
    refs.formError.appendChild(h("span", { class: "ico", "aria-hidden": "true" }, "✖"));
    refs.formError.appendChild(h("span", null, text));
    var detail = (err && (err.userMessage || err.message)) || "";
    refs.formError.appendChild(h("span", { class: "actions" },
      h("button", { type: "button", onclick: function () { run(false); } }, S.errors.tryAgain),
      h("button", { type: "button", onclick: function () { copyText(detail + "\n" + JSON.stringify(state.lastRequests)).catch(function () { announce(S.results.copyFailed); }); } }, S.errors.copyError)));
    showView("form"); announce(text);
  }
  function copyApiRequest() {
    var reqs = L.buildRequests(collectState(), api.defaults);
    copyWithFeedback(JSON.stringify(reqs.length === 1 ? reqs[0] : reqs, null, 2), refs.copyMsg, S.advanced.developerOptions.copied);
  }

  /* =====================================================================
     RESULTS
     ===================================================================== */
  function showView(which) {
    refs.form.hidden = which !== "form"; refs.results.hidden = which !== "results";
    if (which === "results") { var b = refs.results.querySelector(".banner"); if (b) b.focus(); window.scrollTo(0, 0); }
    else window.scrollTo(0, 0);
  }
  function seqTitle(seq) { return hostName(seq.host) + " · " + (S.goalNames[seq.mode] || seq.mode); }

  function renderResults() {
    var report = state.lastReport; clear(refs.results);
    var tabs = null;
    var body = h("div", { id: "result-body" });
    if (report.sequences.length > 1) { body.setAttribute("role", "tabpanel"); body.setAttribute("aria-labelledby", "tab-" + state.activeTab); }
    if (report.sequences.length > 1) {
      tabs = h("div", { class: "tabs", role: "tablist", "aria-label": S.results.summary });
      report.sequences.forEach(function (seq, i) {
        var t = h("button", { type: "button", role: "tab", class: "tab", id: "tab-" + i, "aria-selected": String(i === state.activeTab), "aria-controls": "result-body", tabindex: i === state.activeTab ? "0" : "-1" }, seqTitle(seq));
        t.addEventListener("click", function () { selectTab(i); });
        t.addEventListener("keydown", function (e) {
          var n = report.sequences.length, to = e.key === "ArrowRight" ? (i + 1) % n : e.key === "ArrowLeft" ? (i + n - 1) % n : e.key === "Home" ? 0 : e.key === "End" ? n - 1 : null;
          if (to !== null) { e.preventDefault(); selectTab(to); refs.results.querySelector("#tab-" + to).focus(); }
        });
        tabs.appendChild(t);
      });
    }
    add(refs.results, [tabs, body, report.sequences.length > 1 ? summaryTable(report) : null, resultsFooter()]);
    refs.body = body; renderSequence();
  }
  function selectTab(i) {
    state.activeTab = i;
    Array.prototype.forEach.call(refs.results.querySelectorAll(".tab"), function (t, j) { t.setAttribute("aria-selected", String(j === i)); t.tabIndex = j === i ? 0 : -1; });
    if (refs.body.getAttribute("role") === "tabpanel") refs.body.setAttribute("aria-labelledby", "tab-" + i);
    renderSequence(); announce(seqTitle(state.lastReport.sequences[i]));
  }
  function summaryTable(report) {
    var tb = h("tbody");
    report.sequences.forEach(function (seq, i) {
      var c = L.statusCounts(seq);
      tb.appendChild(h("tr", null,
        h("td", { "data-label": S.results.summaryHost }, hostName(seq.host)),
        h("td", { "data-label": S.results.summaryGoal }, S.goalNames[seq.mode] || seq.mode),
        h("td", { "data-label": S.results.summaryChecks }, statusPill(c.review + c.fail === 0 ? "pass" : c.fail ? "fail" : "review", F(S.results.ofChecks, { p: c.pass, t: c.total })),
          " ", h("button", { type: "button", class: "link small", onclick: function () { selectTab(i); refs.results.querySelector(".banner").focus(); } }, seqTitle(seq)))));
    });
    return h("div", { class: "panel" }, h("h2", null, S.results.summary),
      h("table", { class: "stack" }, h("thead", null, h("tr", null, h("th", null, S.results.summaryHost), h("th", null, S.results.summaryGoal), h("th", null, S.results.summaryChecks))), tb));
  }
  function statusPill(status, text) {
    var ico = status === "pass" ? "✓" : status === "review" ? "⚠" : "✖";
    var word = status === "pass" ? S.results.statusPass : status === "review" ? S.results.statusReview : S.results.statusFail;
    return h("span", { class: "pill " + status }, h("span", { "aria-hidden": "true" }, ico), text ? text : word);
  }
  var CHECK_TERMS = { start_stop: ["stop"], cut_sites: ["cutsite"], rare_codons: ["rare"], codon_score: ["cai"], synthesis: ["synthesis", "repeat", "gc"], start_region: ["rbs"], hairpins: ["hairpin", "fold"] };

  function renderSequence() {
    var report = state.lastReport, seq = report.sequences[state.activeTab], body = clear(refs.body);
    var banner = L.bannerText(seq);
    body.appendChild(h("div", { class: "banner " + (banner.ok ? "ok" : "some"), tabindex: "-1", role: "status" },
      h("span", { "aria-hidden": "true" }, banner.ok ? "✓" : "⚠"), h("span", null, banner.text)));

    // sequence card
    var proteinBox = h("div", { class: "preview", hidden: true, id: "protein-out" });
    var toggle = h("button", { type: "button", "aria-expanded": "false", "aria-controls": "protein-out" }, S.results.showProtein);
    toggle.addEventListener("click", function () {
      var open = proteinBox.hidden; proteinBox.hidden = !open; toggle.setAttribute("aria-expanded", String(open)); toggle.textContent = open ? S.results.hideProtein : S.results.showProtein;
      if (open) proteinBox.textContent = (state.lastRequests && state.lastRequests[0].vector_provides_start ? "M" : "") + L.translateFrame1(seq.dna);
    });
    var copyMsg = h("span", { class: "muted small", role: "status" });
    var tools = h("div", { class: "row seq-tools" },
      h("button", { type: "button", onclick: function () { copyWithFeedback(seq.dna, copyMsg, S.results.copied); } }, S.results.copyDna),
      h("button", { type: "button", onclick: function () { download(L.fastaText(an.name || "optimized_sequence", seq.dna), "optimized_sequence.fasta"); } }, S.results.downloadFasta),
      toggle, copyMsg);
    var dna = h("div", { class: "dna", role: "group", "aria-label": S.results.dnaCaption });
    L.dnaRows(seq.dna, 5).forEach(function (row) {
      dna.appendChild(h("div", { class: "line" }, h("span", { class: "pos", "aria-hidden": "true" }, String(row.start)),
        h("span", { class: "blocks" }, row.blocks.map(function (b) { return h("span", null, b); }))));
    });
    var notes = [];
    if (state.lastRequests && state.lastRequests[0].vector_provides_start) notes.push(S.results.startNote);
    var regionCount = state.lastRequests && state.lastRequests[0].structural_regions ? state.lastRequests[0].structural_regions.length : 0;
    if (regionCount && state.regions.some(function (r) { return r.auto; }) && !state.regionsEdited) notes.push(F(S.results.regionsNote, { n: regionCount }));
    if (state.seedUsed !== null) notes.push(F(S.advanced.seed.used, { seed: state.seedUsed }));
    body.appendChild(h("div", { class: "panel" }, h("h2", null, seqTitle(seq)), tools, dna, proteinBox,
      notes.map(function (n) { return h("p", { class: "help" }, n); })));

    // checklist
    var checkList = function (group) {
      return h("ul", { class: "checklist" }, seq.checks.filter(function (c) { return c.group === group && c.applicable; }).map(function (c) {
        var tpl = S.checks[c.id], det = h("details", null, h("summary", null, S.results.whatThisMeans), h("p", { style: "margin:4px 0" }, tpl.means, " ", h("span", { class: "muted small" }, L.techName(c))));
        var nameEl = h("span", { class: "name" }, tpl.name); (CHECK_TERMS[c.id] || []).slice(0, 1).forEach(function (k) { nameEl.appendChild(tip(k)); });
        return h("li", null, h("span", null, statusPill(c.status)),
          h("div", null, nameEl, h("div", { class: "resultline" }, h("span", null, L.checkResultText(c)), det)));
      }));
    };
    body.appendChild(h("div", { class: "panel" }, h("h2", null, S.results.mustPass), checkList("must_pass"), h("h2", { style: "margin-top:20px" }, S.results.quality), checkList("quality")));

    // fix panel
    var fix = L.fixPanel(seq, report.limits_used);
    if (fix && seq._kept) body.appendChild(msgEl("warn", "\u26A0", S.results.kept));
    else if (fix) body.appendChild(fixPanel(fix, seq));
    body.appendChild(partFinderPanel(seq));
  }
  function partFinderPanel(seq) {
    var P = S.parts, out = h("div", { class: "parts-out" }), status = h("p", { class: "muted", role: "status" });
    var sel = h("select", { id: "pf-level" });
    ["any", "low", "moderate", "high"].forEach(function (k) { sel.appendChild(h("option", { value: k }, P.levels[k])); });
    var btn = h("button", { type: "button", class: "primary" }, P.find);
    btn.addEventListener("click", function () {
      if (btn.getAttribute("aria-disabled") === "true") return;
      btn.setAttribute("aria-disabled", "true"); status.textContent = P.finding; clear(out);
      var body = { host: seq.host, top: 5, cds: seq.dna }; if (sel.value !== "any") body.level = sel.value;
      fetch("/api/part-finder", { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(body) })
        .then(function (r) { return r.json().then(function (d) { return { ok: r.ok, d: d }; }); })
        .then(function (x) {
          btn.removeAttribute("aria-disabled"); status.textContent = "";
          if (!x.ok) { out.appendChild(msgEl("error", "\u2716", /chassis tag/.test(x.d.error || "") ? P.unsupported : (x.d.error || P.error))); return; }
          renderParts(out, x.d, seq);
        })
        .catch(function () { btn.removeAttribute("aria-disabled"); status.textContent = ""; out.appendChild(msgEl("error", "\u2716", P.error)); });
    });
    return h("div", { class: "panel", id: "part-finder" }, h("h2", null, P.title), h("p", { class: "help" }, P.intro),
      msgEl("info", "\u2139", P.honest),
      h("div", { class: "field" }, h("label", { class: "label", for: "pf-level" }, P.levelLabel), h("p", { class: "help" }, P.levelHelp), sel),
      h("div", { class: "row" }, btn, status), out);
  }
  function buildPanel(partsOut, seq) {
    var B = S.build, result = h("div", { class: "build-result" }), status = h("p", { class: "muted", role: "status" });
    var btn = h("button", { type: "button", class: "primary" }, B.button);
    btn.addEventListener("click", function () {
      if (btn.getAttribute("aria-disabled") === "true") return;
      var pick = function (k) { var c = partsOut.querySelector("input[name=pf-" + k + "]:checked"); return c ? c.value : null; };
      var promoter = pick("promoter"), rbs = pick("rbs"), terminator = pick("terminator");
      clear(result);
      if (!promoter || !rbs || !terminator) { result.appendChild(msgEl("warn", "\u26A0", B.needParts)); return; }
      btn.setAttribute("aria-disabled", "true"); status.textContent = B.building;
      var vectorAtg = state.lastRequests && state.lastRequests[0].vector_provides_start;
      fetch("/api/vectorize", { method: "POST", headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ promoter: promoter, rbs: rbs, terminator: terminator, cds: (vectorAtg ? "ATG" : "") + seq.dna }) })
        .then(function (r) { return r.json().then(function (d) { return { ok: r.ok, d: d }; }); })
        .then(function (x) {
          btn.removeAttribute("aria-disabled"); status.textContent = "";
          if (!x.ok) { result.appendChild(msgEl("error", "\u2716", x.d.error || B.error)); return; }
          renderBuild(result, x.d);
        })
        .catch(function () { btn.removeAttribute("aria-disabled"); status.textContent = ""; result.appendChild(msgEl("error", "\u2716", B.error)); });
    });
    return h("div", { class: "field", style: "border-top:1px solid var(--line);padding-top:16px;margin-top:24px" },
      h("h2", null, B.title), h("p", { class: "help" }, B.intro), h("div", { class: "row" }, btn, status), result);
  }
  function renderBuild(out, d) {
    var B = S.build, g = d.gc, s = d.size, j = d.junction_check;
    out.appendChild(h("p", { class: "help" }, F(B.backbone, { name: d.backbone.name, bp: d.backbone.length, desc: d.backbone.description })));
    if (j.ok) out.appendChild(msgEl("note", "\u2713", B.junctionOk));
    else {
      out.appendChild(msgEl("error", "\u2716", B.junctionBad));
      out.appendChild(h("ul", null, j.violations.map(function (v) { return h("li", null, F(B.violation, { enzyme: v.enzyme, start: v.start, end: v.end, where: v.at_junction ? B.atJoin : "" })); })
        .concat(j.type_iis_sites.map(function (n) { return h("li", null, n); }))));
    }
    out.appendChild(h("p", { class: "help" }, B.flankNote));
    out.appendChild(h("h3", { style: "margin-top:16px" }, B.sizeTitle));
    out.appendChild(h("p", null, F(B.size, { total: s.total_bp.toLocaleString("en-US"), backbone: s.backbone_bp.toLocaleString("en-US"), insert: s.insert_bp.toLocaleString("en-US") }), " ", s.within_studied_range ? B.sizeWithin : B.sizeBeyond));
    out.appendChild(h("p", { class: "help" }, s.evidence));
    if (d.synthesis.beyond_clonal_gene_limit) out.appendChild(msgEl("warn", "\u26A0", F(B.synthBad, { limit: "7,000" })));
    out.appendChild(h("h3", { style: "margin-top:16px" }, B.gcTitle));
    out.appendChild(h("p", null, F(B.gc, { gc: g.insert_percent, lo: g.window_min, hi: g.window_max })));
    out.appendChild(g.twist_high_complexity ? msgEl("warn", "\u26A0", B.gcTwistBad) : msgEl("note", "\u2713", B.gcTwistOk));
    if (g.outside_project_range) out.appendChild(msgEl("info", "\u2139", B.gcProject));
    out.appendChild(h("h3", { style: "margin-top:16px" }, B.layoutTitle));
    out.appendChild(h("table", { class: "stack" }, h("tbody", null, d.layout.map(function (x) {
      return h("tr", null, h("td", { "data-label": "Piece" }, B.kind[x.kind]), h("td", { "data-label": "Name" }, x.kind === "scar" ? d.insert.substr(x.start - d.backbone.length - 1, x.length) : x.name),
        h("td", { "data-label": "Positions" }, x.start + "-" + x.end), h("td", { "data-label": "Length" }, String(x.length)));
    }))));
    out.appendChild(h("p", { class: "help" }, F(B.scarNote, { std: d.scars.standard, rbs: d.scars.rbs_to_cds })));
    var msg = h("span", { class: "muted small", role: "status" });
    out.appendChild(h("h3", { style: "margin-top:16px" }, B.seqTitle));
    out.appendChild(h("div", { class: "row seq-tools" },
      h("button", { type: "button", onclick: function () { copyWithFeedback(d.sequence, msg, B.copied); } }, B.copy),
      h("button", { type: "button", onclick: function () { download(d.fasta, "plasmid_" + d.backbone.name + ".fasta"); } }, B.download), msg));
    var dna = h("div", { class: "dna", role: "group", "aria-label": B.seqTitle });
    L.dnaRows(d.sequence, 5).forEach(function (row) {
      dna.appendChild(h("div", { class: "line" }, h("span", { class: "pos", "aria-hidden": "true" }, String(row.start)), h("span", { class: "blocks" }, row.blocks.map(function (b) { return h("span", null, b); }))));
    });
    out.appendChild(h("details", null, h("summary", null, B.seqTitle + " (" + d.sequence.length.toLocaleString("en-US") + ")"), dna));
    out.appendChild(h("p", { class: "help" }, B.sources + ": ", d.sources.scars.concat([d.sources.synthesis]).map(function (u, i) { return h("a", { href: u, target: "_blank", rel: "noopener", style: "margin-right:12px" }, B.sourceLabels[i]); })));
  }
  function renderParts(out, d, seq) {
    var P = S.parts;
    if (d.cds_rfc10) out.appendChild(d.cds_rfc10.ok ? msgEl("note", "\u2713", P.cdsOk) : msgEl("error", "\u2716", F(P.cdsBad, { sites: d.cds_rfc10.illegal_sites.join(", ") })));
    (d.warnings || []).forEach(function (w) { out.appendChild(msgEl("warn", "\u26A0", w)); });
    ["promoter", "rbs", "terminator"].forEach(function (kind) {
      var block = d.results[kind]; if (!block) return;
      out.appendChild(h("h3", { style: "margin-top:20px" }, P.kinds[kind]));
      out.appendChild(h("p", { class: "help" }, F(P.counts, { p: block.passing_filters, c: block.considered })));
      if (!block.parts.length) { out.appendChild(h("p", { class: "muted" }, P.none)); return; }
      out.appendChild(h("ul", { class: "partlist" }, block.parts.map(function (part) {
        var avail = part.status === "Available", e = part.expression;
        var strength = e.verified_strength ? F(P.strengthKnown, { value: e.verified_strength.value, unit: e.verified_strength.unit }) : P.strengthNone;
        var lines = [strength]; if (part.regulation) lines.unshift(P.regulation[part.regulation]);
        if (e.level_hint) lines.push(F(P.hintText, { word: "\u2018" + e.level_hint.word + "\u2019" }));
        return h("li", { class: "part" },
          h("label", { class: "check", style: "min-height:44px;padding:0" },
            h("input", { type: "radio", name: "pf-" + kind, value: part.name, checked: block.parts[0] === part }), S.build.use),
          h("div", { class: "row between", style: "margin:0" },
            h("a", { href: part.registry_url, target: "_blank", rel: "noopener" }, part.name, h("span", { class: "sr-only" }, " " + P.viewRegistry)),
            h("span", { class: "pill " + (avail ? "pass" : "review") }, h("span", { "aria-hidden": "true" }, avail ? "\u2713" : "\u26A0"), avail ? P.available : P.notAvailable)),
          part.description ? h("p", { class: "help", style: "margin:2px 0" }, part.description) : null,
          h("p", { class: "small", style: "margin:2px 0" }, part.chassis_evidence === "tagged" ? P.hostTagged : P.hostUnknown, " \u00B7 ", lines.join(" \u00B7 ")),
          part.sequence.length <= 60 ? h("p", { class: "dna small", style: "margin:2px 0;overflow-wrap:anywhere" }, part.sequence) : null,
          h("details", null, h("summary", null, P.why), h("ul", null, part.reasons.map(function (r) { return h("li", null, r); }))));
      })));
    });
    out.appendChild(buildPanel(out, seq));
    out.appendChild(h("p", { class: "help", style: "margin-top:16px" }, P.limits));
    out.appendChild(h("p", { class: "help" }, F(P.source, { file: d.dataset.source.file, sha: d.dataset.source.sha256.slice(0, 12) })));
  }
  function fixPanel(fix, seq) {
    var msg = h("p", { class: "help", role: "status" });
    var row = h("div", { class: "row" }, fix.buttons.map(function (b) {
      return h("button", { type: "button", onclick: function () {
        if (b.keep) { seq._kept = true; renderSequence(); announce(S.results.kept); return; }
        state.overrides = Object.assign({}, state.overrides, b.override); showView("form"); run(false);
      } }, b.label);
    }));
    return h("div", { class: "fix", role: "region", "aria-label": S.results.fixTitle },
      h("h3", null, "⚠ " + S.results.fixTitle), h("blockquote", null, fix.text), row, msg);
  }
  function resultsFooter() {
    var report = state.lastReport;
    var tech = h("dl", { class: "techlist" });
    var seq = function () { return report.sequences[state.activeTab]; };
    var techBox = h("details", { class: "panel" }, h("summary", { class: "label" }, S.results.fullReport), tech);
    techBox.addEventListener("toggle", function () {
      if (!techBox.open) return; clear(tech); var s = seq(), T = S.tech;
      function row(k, v) { tech.appendChild(h("dt", null, k)); tech.appendChild(h("dd", null, String(v))); }
      row(T.cai, s.cai === null ? T.none : s.cai.toFixed(3)); row(T.tai, s.tai === null ? T.taiNone : s.tai.toFixed(3));
      row(T.minW, s.min_w_used.toFixed(3)); row(T.belowW, s.codons_below_w_threshold);
      row(T.startDg, s.five_prime_dG === null ? T.none : s.five_prime_dG.toFixed(2) + " " + T.kcal);
      row(T.startOpen, s.init_region_unpaired === null ? T.none : (s.init_region_unpaired ? T.yes : T.no));
      row(T.worstMfe, s.worst_window_mfe === null ? T.none : s.worst_window_mfe.toFixed(2) + " " + T.kcal);
      row(T.stem, s.longest_stem); row(T.terminators, s.terminator_motifs.length); row(T.inverted, s.inverted_repeats.length);
      row(T.gc, (s.gc_overall * 100).toFixed(1) + "%"); row(T.gcWindow, (s.gc_window_min * 100).toFixed(0) + "% – " + (s.gc_window_max * 100).toFixed(0) + "%");
      row(T.repeats, s.repeated_15mers); row(T.pauses, s.pause_site_count);
      var src = report.codon_table_sources[s.host]; if (src) row(T.sources, src.source + " (" + src.retrieved + ")");
      s.notes.forEach(function (n) { row(T.notes, n); });
    });
    return h("footer", { class: "page" },
      h("div", { class: "row" },
        h("button", { type: "button", onclick: function () { showView("form"); refs.ta.focus(); } }, S.results.editSettings),
        h("button", { type: "button", onclick: function () { showView("form"); run(false); } }, S.results.runAgain)),
      techBox);
  }
  function download(text, name) {
    var url = URL.createObjectURL(new Blob([text], { type: "text/plain" }));
    var a = h("a", { href: url, download: name }); document.body.appendChild(a); a.click(); document.body.removeChild(a); setTimeout(function () { URL.revokeObjectURL(url); }, 1000);
  }

  /* =====================================================================
     BOOT
     ===================================================================== */
  function boot() {
    document.title = S.app.title;
    var app = document.getElementById("app");
    refs.live = h("div", { class: "sr-only", "aria-live": "polite", "aria-atomic": "true", id: "live" });
    refs.results = h("div", { id: "resultsView", hidden: true });
    var seen = false; try { seen = localStorage.getItem("pd_seen") === "1"; } catch (e) { /* ignore */ }
    add(app, [h("header", { class: "top" }, h("h1", null, S.app.title), h("a", { href: "/glossary" }, S.app.glossaryLink)),
      buildForm(seen), refs.results, h("footer", { class: "page" }, S.app.privacy), refs.live]);
    if (load("pd_adv") === "1") setAdvanced(true);
    onProteinChange(); onGoalChange(); renderChips(); updateBadge(); showLimitDefaults(); refreshStartParts();
    fetch("/api/limits").then(function (r) { return r.json(); }).then(function (d) { api = d; onProteinChange(); syncTemperature(); showLimitDefaults(); }).catch(function () { /* keep built-in defaults */ });
  }
  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", boot); else boot();
})();
