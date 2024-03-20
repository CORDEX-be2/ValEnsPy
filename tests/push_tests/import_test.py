import sys, os
from pathlib import Path

lib_folder = Path(__file__).resolve().parents[2].joinpath("src")

import valenspy

print(f"Successfully imported valenspy version: {valenspy.__version__}")