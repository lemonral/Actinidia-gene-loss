from __future__ import annotations

import importlib.util
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "figures" / "render_manuscript_unit_loss.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("render_manuscript_unit_loss", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_publication_labels_use_counts_and_grey_deleted_calls() -> None:
    module = _load_module()
    assert module.PUBLICATION_VALUE_COLUMN == "gene_count"
    assert module.PUBLICATION_CATEGORY_COLORS["shared_deleted"] == "#595959"
    assert module.PUBLICATION_CATEGORY_COLORS["nonshared_deleted"] == "#B3B3B3"
    assert module.PUBLICATION_CATEGORY_LABELS == {
        "shared_decayed": "Shared decayed",
        "shared_deleted": "Shared deleted",
        "nonshared_decayed": "Non-shared decayed",
        "nonshared_deleted": "Non-shared deleted",
    }
    visible_text = " ".join(module.PUBLICATION_CATEGORY_LABELS.values()).lower()
    assert "article" not in visible_text
    assert "threshold" not in visible_text
