# script to extract onsite density matrices and use them to construct h5 archive with vasp dft_tools interface for nickel 2 band eg model

# load needed modules
import numpy as np

import os.path
import shutil
import re
import cmath

from triqs_dft_tools.sumk_dft import *
from triqs_dft_tools.sumk_dft_tools import *
from triqs_dft_tools.converters.vasp import *
from triqs_dft_tools.converters.plovasp.vaspio import VaspData
from triqs_dft_tools.converters.plovasp.plotools import generate_plo, output_as_text
import triqs_dft_tools.converters.plovasp.converter as plo_converter


vasp_dir = './'

# Generate and store PLOs
plo_converter.generate_and_output_as_text('plo.cfg', vasp_dir='./')

# run the archive creat routine
conv = VaspConverter('vasp')
conv.convert_dft_input()

sk = SumkDFTTools(hdf_file='vasp.h5')
