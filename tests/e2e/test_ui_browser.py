"""Real-browser checks for the redesigned UI. Skipped unless Playwright and a
Chromium build are installed (pip install playwright; playwright install chromium).
Runs the Flask app in-process on a free port. axe-core is loaded from a local
file named by AXE_PATH or from cdnjs; the accessibility test skips if neither works."""
import os
import socket
import threading
import urllib.request

import pytest

sync_api = pytest.importorskip("playwright.sync_api")

from api.index import app  # noqa: E402

UBQ = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
AXE_CDN = "https://cdnjs.cloudflare.com/ajax/libs/axe-core/4.10.2/axe.min.js"


@pytest.fixture(scope="module")
def base_url():
    s = socket.socket(); s.bind(("127.0.0.1", 0)); port = s.getsockname()[1]; s.close()
    from werkzeug.serving import make_server
    server = make_server("127.0.0.1", port, app, threaded=True)
    threading.Thread(target=server.serve_forever, daemon=True).start()
    yield f"http://127.0.0.1:{port}"
    server.shutdown()


@pytest.fixture(scope="module")
def browser():
    with sync_api.sync_playwright() as p:
        try:
            b = p.chromium.launch()
        except Exception as exc:  # browser not installed
            pytest.skip(f"Chromium unavailable: {exc}")
        yield b
        b.close()


@pytest.fixture()
def page(browser, base_url):
    ctx = browser.new_context(viewport={"width": 1100, "height": 900})
    pg = ctx.new_page()
    pg.errors = []
    pg.on("pageerror", lambda e: pg.errors.append(str(e)))
    pg.on("console", lambda m: pg.errors.append(m.text) if m.type == "error" else None)
    pg.goto(base_url + "/")
    pg.wait_for_selector("#protein")
    yield pg
    assert pg.errors == []
    ctx.close()


def paste(page, text):
    page.fill("#protein", text)


def test_default_path_is_three_clicks_or_fewer(page):
    """Acceptance #1: after pasting, a result takes at most 3 clicks."""
    paste(page, UBQ)
    page.click("#optimize")
    page.wait_for_selector(".banner", timeout=120000)
    assert page.locator(".checklist li").count() == 10


def test_optimize_disabled_until_valid(page):
    assert page.get_attribute("#optimize", "aria-disabled") == "true"
    page.click("#optimize", force=True)  # aria-disabled buttons stay clickable so they can explain themselves
    assert "Paste a protein sequence to start." in page.inner_text(".actions-bar")
    assert page.locator(".banner").count() == 0
    paste(page, UBQ)
    assert page.get_attribute("#optimize", "aria-disabled") == "false"


def test_invalid_letters_are_highlighted_and_removable(page):
    paste(page, "MKTAY-IAKQRQ?ISFVKSHFSRQ")
    assert "Found 2 letters that aren't amino acids (highlighted)" in page.inner_text("#protein-msgs")
    assert page.locator(".preview mark").count() == 2
    page.click("text=Remove them")
    assert page.input_value("#protein") == "MKTAYIAKQRQISFVKSHFSRQ"


def test_dna_paste_offers_translation(page):
    paste(page, "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
    page.click("text=Translate it for me")
    assert page.input_value("#protein") == "MAIVMGR"


def test_no_host_message_and_focus(page):
    paste(page, UBQ)
    page.uncheck("input[value=e_coli_bl21_de3]")
    page.click("#optimize")
    assert "Pick at least one place to make your protein." in page.inner_text("fieldset")


def test_host_search_and_more_organisms(page):
    page.click("summary:has-text('More organisms')")
    page.fill("#host-search", "bacillus")
    visible = page.locator(".card:has(input[type=checkbox]):visible").count()
    assert visible == 2
    page.fill("#host-search", "zzz")
    assert page.locator("text=No organism matches that search.").is_visible()


def test_internal_ids_never_shown_on_default_view(page):
    text = page.inner_text("#app")
    for hid in ("e_coli_bl21_de3", "s_cerevisiae", "b_subtilis_168", "{", "}"):
        assert hid not in text


def test_advanced_badge_counts_changes_and_resets(page):
    page.click("text=Advanced settings")
    page.check("#hedge")
    page.fill("#seed", "42")
    assert "Advanced settings (2 changed)" in page.inner_text("#app")
    page.click("text=Reset to defaults")
    assert "changed" not in page.inner_text(".adv-toggle")


def test_copy_api_request_sends_only_changed_values(page):
    page.context.grant_permissions(["clipboard-read", "clipboard-write"])
    paste(page, UBQ)
    page.click("text=Advanced settings")
    page.fill("#seed", "42")
    page.click("summary:has-text('Developer options')")
    page.click("text=Copy API request")
    import json
    body = json.loads(page.evaluate("navigator.clipboard.readText()"))
    assert body == {"protein": UBQ, "mode": "BALANCED", "hosts": ["e_coli_bl21_de3"], "seed": 42}


def test_multi_host_gives_tabs_and_summary(page):
    paste(page, UBQ)
    page.check("input[value=human]")
    page.click("#optimize")
    page.wait_for_selector(".tab", timeout=180000)
    assert page.locator(".tab").count() == 2
    page.click(".tab >> nth=1")
    assert "Human cells" in page.inner_text("#result-body h2 >> nth=0")
    assert page.locator("table.stack tbody tr").count() == 2


def test_fix_panel_offers_explicit_buttons_and_reuses_seed(page):
    """A hairpin the engine can't fully remove must surface as a fix panel
    with buttons, and re-running must use the same number as before."""
    paste(page, UBQ)
    page.click("#optimize")
    page.wait_for_selector(".banner", timeout=120000)
    seed_text = page.inner_text("text=Number used")
    if page.locator(".fix").count() == 0:
        pytest.skip("this run happened to pass every check")
    assert "Keep this result" in page.inner_text(".fix")
    page.click(".fix button:has-text('re-run')")
    page.wait_for_selector(".banner", timeout=120000)
    assert page.inner_text("text=Number used") == seed_text


def test_escape_closes_tooltip(page):
    page.click("legend .tip")
    assert page.locator("[role=tooltip]:visible").count() == 1
    page.keyboard.press("Escape")
    assert page.locator("[role=tooltip]:visible").count() == 0


def test_keyboard_only_default_flow(page):
    page.focus("#protein")
    page.keyboard.type(UBQ)
    for _ in range(40):  # tab to the Optimize button
        page.keyboard.press("Tab")
        if page.evaluate("document.activeElement.id") == "optimize":
            break
    assert page.evaluate("document.activeElement.id") == "optimize"
    page.keyboard.press("Enter")
    page.wait_for_selector(".banner", timeout=120000)


def test_goal_cards_are_a_radio_group_with_arrow_keys(page):
    page.focus("input[name=goal]:checked")
    page.keyboard.press("ArrowRight")
    assert page.evaluate("document.querySelector('input[name=goal]:checked').value") != "BALANCED"


def test_phone_layout_has_sticky_button_and_no_horizontal_scroll(browser, base_url):
    ctx = browser.new_context(viewport={"width": 390, "height": 800}, is_mobile=True)
    pg = ctx.new_page(); pg.goto(base_url + "/"); pg.wait_for_selector("#protein")
    assert pg.evaluate("document.documentElement.scrollWidth <= window.innerWidth")
    assert pg.evaluate("getComputedStyle(document.querySelector('.actions-bar')).position") == "sticky"
    cols = pg.evaluate("getComputedStyle(document.querySelector('.cards.four')).gridTemplateColumns.split(' ').length")
    assert cols == 1
    ctx.close()


def test_touch_targets_at_least_44px(page):
    small = page.evaluate("""() => [...document.querySelectorAll('button, input[type=checkbox], input[type=radio], select, summary')]
        .filter(e => e.offsetParent !== null && !e.closest('.tipwrap') && !e.classList.contains('tip'))
        .map(e => (e.matches('input') && e.closest('label') ? e.closest('label') : e))  // the label is the tap target
        .map(e => { const r = e.getBoundingClientRect(); return [e.tagName + ' ' + (e.textContent || e.name || '').slice(0, 20), Math.round(r.width), Math.round(r.height)]; })
        .filter(([, w, h]) => h < 44 || w < 44)""")
    assert small == []


@pytest.mark.parametrize("scheme", ["light", "dark"])
def test_axe_no_critical_or_serious_issues(browser, base_url, scheme):
    path = os.environ.get("AXE_PATH")
    try:
        axe = open(path).read() if path and os.path.exists(path) else urllib.request.urlopen(AXE_CDN, timeout=20).read().decode()
    except Exception as exc:
        pytest.skip(f"axe-core unavailable: {exc}")
    ctx = browser.new_context(viewport={"width": 1100, "height": 900}, color_scheme=scheme)
    pg = ctx.new_page(); pg.goto(base_url + "/"); pg.wait_for_selector("#protein")

    def scan(label):
        pg.evaluate(axe)
        res = pg.evaluate("axe.run(document, {runOnly: {type: 'tag', values: ['wcag2a','wcag2aa','wcag21a','wcag21aa','wcag22aa']}})")
        bad = [(v["id"], v["impact"], [n["target"] for n in v["nodes"]][:3]) for v in res["violations"] if v["impact"] in ("critical", "serious")]
        assert bad == [], f"{label} ({scheme}): {bad}"

    scan("form")
    pg.click("text=Try an example"); pg.click("text=Advanced settings"); pg.click("summary:has-text('Developer options')")
    scan("form with advanced open")
    pg.click("#optimize"); pg.wait_for_selector(".banner", timeout=120000)
    scan("results")
    pg.click("#part-finder button.primary"); pg.wait_for_selector(".partlist", timeout=30000)
    scan("results with part finder")
    ctx.close()


# ---- Part Finder and the button fixes ----------------------------------------------------------
def run_example(page):
    paste(page, UBQ)
    page.click("#optimize")
    page.wait_for_selector(".banner", timeout=120000)


def test_start_of_gene_options_come_from_the_registry_extract(page):
    page.click("text=Advanced settings")
    page.wait_for_function("document.querySelector('#start-choice optgroup').children.length > 3")
    opts = page.locator("#start-choice optgroup option").all_inner_texts()
    assert "BBa_B0034 \u2014 AAAGAGGAGAAA" in opts
    page.select_option("#start-choice", "part:BBa_B0034")
    page.context.grant_permissions(["clipboard-read", "clipboard-write"])
    paste(page, UBQ)
    page.click("summary:has-text('Developer options')")
    page.click("text=Copy API request")
    import json
    assert json.loads(page.evaluate("navigator.clipboard.readText()"))["five_prime_utr"] == "AAAGAGGAGAAA"


def test_part_finder_lists_three_kinds_with_honest_notes_and_provenance(page):
    run_example(page)
    page.click("#part-finder button.primary")
    page.wait_for_selector(".partlist", timeout=30000)
    panel = page.inner_text("#part-finder")
    for heading in ("Promoters", "Ribosome binding sites", "Terminators"):
        assert heading in panel
    assert "can't predict how much protein a part makes" in panel
    assert "No measured strength on file" in panel
    assert "Always on (no inducer needed)" in panel
    assert "checksum c64bbc9a1" in panel
    assert "Your gene passes the RFC10 check." in panel
    hrefs = page.locator("#part-finder .partlist a").evaluate_all("els => els.map(e => e.href)")
    assert hrefs and all(h.startswith("https://registry.igem.org/parts/bba-") for h in hrefs)
    assert page.locator("#part-finder .partlist > li.part").count() <= 15


def test_part_finder_explains_why_a_part_was_chosen(page):
    run_example(page)
    page.click("#part-finder button.primary"); page.wait_for_selector(".partlist")
    page.locator("#part-finder details >> nth=0").locator("summary").click()
    assert "Passes the RFC10 check" in page.locator("#part-finder details >> nth=0").inner_text()


def test_part_finder_says_plainly_when_a_host_has_no_registry_tag(page):
    page.uncheck("input[value=e_coli_bl21_de3]")
    page.click("summary:has-text('More organisms')")
    page.check("input[value=c_reinhardtii]")
    run_example(page)
    page.click("#part-finder button.primary")
    page.wait_for_selector("#part-finder .msg.error")
    assert "don't have Registry host tags for this organism" in page.inner_text("#part-finder")
    assert any("400" in e for e in page.errors)       # the browser logs the deliberate 400; nothing else may be logged
    page.errors[:] = [e for e in page.errors if "400" not in e]


def test_keep_this_result_shows_a_visible_confirmation(page):
    paste(page, UBQ)
    page.click("text=Advanced settings")
    page.fill("#lim-mfe", "-3")                       # an impossible fold limit forces a conflict
    page.click("#optimize")
    page.wait_for_selector(".fix", timeout=120000)
    page.click(".fix >> text=Keep this result")
    assert page.locator(".fix").count() == 0
    assert "Kept. This result still has the problem described above." in page.inner_text("#result-body")


def test_stale_error_clears_once_the_problem_is_fixed(page):
    page.click("#optimize", force=True)
    assert page.locator("#optimize-note ~ .msg.error").is_visible()
    paste(page, UBQ)
    assert page.locator("#optimize-note ~ .msg.error").is_hidden()


def test_copy_failure_is_reported_not_swallowed(page):
    run_example(page)
    page.evaluate("Object.defineProperty(navigator, 'clipboard', {value: {writeText: () => Promise.reject(new Error('denied'))}, configurable: true})")
    page.click("text=Copy DNA")
    page.wait_for_selector("text=We couldn't copy that automatically")


def test_remove_region_buttons_have_distinct_names(page):
    paste(page, "MKTAYIAKQRQISFVKS" + "GGGGSGGGGS" + "HFSRQLEERLGLIEVQ" + "SSSSSGGGGG" + "APILSRVGDGTQ")
    page.click("text=Advanced settings")
    labels = page.locator("table.regions button").evaluate_all("els => els.map(e => e.getAttribute('aria-label'))")
    assert labels == ["Remove region 1", "Remove region 2"]


def test_loading_state_takes_focus_and_can_be_cancelled(page):
    page.route("**/api/optimize", lambda route: None)   # hold the request open so the loading state is stable
    paste(page, UBQ)
    page.click("#optimize")
    page.wait_for_selector(".loading:not([hidden])")
    assert page.evaluate("document.activeElement.classList.contains('loading')")
    page.click(".loading button:has-text('Cancel')")
    assert page.evaluate("document.activeElement.id") == "optimize"


def test_result_tabs_have_a_labelled_tab_panel(page):
    paste(page, UBQ)
    page.check("input[value=human]")
    page.click("#optimize"); page.wait_for_selector(".tab", timeout=180000)
    assert page.get_attribute("#result-body", "role") == "tabpanel"
    assert page.get_attribute("#result-body", "aria-labelledby") == "tab-0"
    page.click(".tab >> nth=1")
    assert page.get_attribute("#result-body", "aria-labelledby") == "tab-1"
