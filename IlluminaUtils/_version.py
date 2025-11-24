from pathlib import Path

def _read_version():
    here = Path(__file__).resolve()
    version_file = here.with_name("VERSION")
    return version_file.read_text(encoding="utf-8").strip()

__version__ = _read_version()
