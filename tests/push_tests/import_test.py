import sys, os
from pathlib import Path

# ----- Import test as package -----#

#This is to test if the package is in your PYTHONPATH, and can be imported
import valenspy

# ----- Import valenspy (NOT AS A PACKAGE) ---------
lib_folder = Path(__file__).resolve().parents[2].joinpath("src")
sys.path.insert(0, str(lib_folder))
import valenspy

print(f"Successfully imported valenspy version: {valenspy.__version__}")
