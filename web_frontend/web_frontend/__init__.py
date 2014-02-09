import logging
import os
import sys


# Append EbeCollector module path to PYTHONPATH.
lib_path = os.path.abspath("../EBE-Node/EbeCollector")
if lib_path not in sys.path:
    sys.path.append(lib_path)

# Log to stdout.
logging.basicConfig(file=sys.stdout, level=logging.INFO)
