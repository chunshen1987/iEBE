import os
import sys

lib_path = os.path.abspath("../EBE-Node/EbeCollector")
if lib_path not in sys.path:
    sys.path.append(lib_path)

