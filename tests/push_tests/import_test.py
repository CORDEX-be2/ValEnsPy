import sys, os
from pathlib import Path

lib_folder = Path(__file__).resolve().parents[2].joinpath("src")

sys.path.insert(0, str(lib_folder))
import valenspy

print(f"Successfully imported valenspy version: {valenspy.__version__}")