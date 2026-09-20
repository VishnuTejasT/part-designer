"""Runs the node unit tests for the browser client (static/logic.js, strings.js).
Uses node's built-in runner, so no npm install is needed."""
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.skipif(shutil.which("node") is None, reason="node is not installed")
def test_client_logic_and_copy_checks():
    r = subprocess.run(["node", "--test", str(ROOT / "tests" / "js" / "logic.test.js")],
                       capture_output=True, text=True, cwd=ROOT, timeout=120)
    assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-1000:]
