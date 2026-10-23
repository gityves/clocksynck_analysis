"""Shares the loaders and local helpers of cancer/_shared.py (same module, loaded under another name)."""

import importlib.util
from pathlib import Path

from common.paths import group_dir  # noqa: E402
_f = group_dir("cancer") / "_shared.py"
_spec = importlib.util.spec_from_file_location("cancer_shared", _f)
_m = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_m)
globals().update({k: v for k, v in vars(_m).items() if not k.startswith("__")})
