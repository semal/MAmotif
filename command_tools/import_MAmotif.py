import os
import sys

current_script_folder = os.path.split(os.path.abspath(__file__))[0]
one_back_folder = os.sep.join(current_script_folder.split(os.sep)[:-1])

sys.path.append(one_back_folder)
from MAmotif_fisher_test import *
from MAmotif_pkg import *