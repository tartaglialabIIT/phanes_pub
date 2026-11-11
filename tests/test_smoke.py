"""Smoke tests that do not import the Dash app (requires profet *.pkl data files)."""

from pathlib import Path


def test_requirements_file_exists():
    req = Path(__file__).resolve().parent.parent / "dash_pip_requirements.txt"
    assert req.is_file()
    text = req.read_text(encoding="utf-8")
    assert "dash" in text.lower()


def test_data_module_loads_without_pickles():
    """data.py only defines loaders; importing it must not touch disk."""
    import data  # noqa: F401

    assert callable(data.generate_scalar_df)
