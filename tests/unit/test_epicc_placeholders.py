"""Unit tests for the `epicc` profile placeholder check.

The placeholder safety net (`check_profile_placeholders`) must catch a
profile that still carries an *active* `<your-...>` value, but must NOT
fire on a placeholder that lives in a commented-out example line — most
notably the shipped `#- slurm_account="<your-slurm-account>"`, which users
are meant to leave alone when their cluster assigns a default SLURM
account. The `epicc` script has no .py extension, so we load it via
importlib (same pattern as test_epicc_output.py).
"""

import importlib.util
import sys
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parents[2]
_EPICC_SRC = _REPO_ROOT / "epicc"


def _load_epicc():
    from importlib.machinery import SourceFileLoader
    loader = SourceFileLoader("epicc_cli", str(_EPICC_SRC))
    spec = importlib.util.spec_from_loader("epicc_cli", loader)
    mod = importlib.util.module_from_spec(spec)
    sys.modules.setdefault("epicc_cli", mod)
    loader.exec_module(mod)
    return mod


epicc = _load_epicc()


# ---------------------------------------------------------------------------
# _strip_yaml_comment
# ---------------------------------------------------------------------------

class TestStripYamlComment:
    def test_full_line_comment_collapses(self):
        line = '  #- slurm_account="<your-slurm-account>"   # EDIT'
        assert epicc._PLACEHOLDER_RE.findall(epicc._strip_yaml_comment(line)) == []

    def test_active_value_with_trailing_comment_preserved(self):
        line = '  - slurm_account="<your-slurm-account>"   # EDIT'
        stripped = epicc._strip_yaml_comment(line)
        assert epicc._PLACEHOLDER_RE.findall(stripped) == ["<your-slurm-account>"]

    def test_hash_inside_quotes_not_treated_as_comment(self):
        line = '  - foo="a#b<your-thing>"'
        assert epicc._strip_yaml_comment(line) == line

    def test_line_without_comment_unchanged(self):
        line = "  - runtime=\"1h\""
        assert epicc._strip_yaml_comment(line) == line


# ---------------------------------------------------------------------------
# check_profile_placeholders
# ---------------------------------------------------------------------------

class TestCheckProfilePlaceholders:
    def _write(self, tmp_path, body):
        cfg = tmp_path / "config.yaml"
        cfg.write_text(body)
        return cfg

    def test_commented_placeholder_accepted(self, tmp_path):
        cfg = self._write(
            tmp_path,
            'default-resources:\n'
            '  - runtime="1h"\n'
            '  #- slurm_account="<your-slurm-account>"   # EDIT\n',
        )
        # Must not raise SystemExit.
        epicc.check_profile_placeholders(str(cfg))

    def test_active_placeholder_rejected(self, tmp_path):
        cfg = self._write(
            tmp_path,
            'default-resources:\n'
            '  - slurm_account="<your-slurm-account>"\n',
        )
        with pytest.raises(SystemExit) as exc:
            epicc.check_profile_placeholders(str(cfg))
        assert exc.value.code == 2

    def test_shipped_slurm_profile_accepted(self):
        # The real shipped example profile leaves slurm_account commented.
        epicc.check_profile_placeholders(str(_REPO_ROOT / "profiles" / "slurm"))

    def test_missing_profile_is_noop(self, tmp_path):
        epicc.check_profile_placeholders(str(tmp_path / "does_not_exist.yaml"))
