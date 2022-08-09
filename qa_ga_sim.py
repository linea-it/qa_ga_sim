import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import Table
import astropy.io.fits as fits
from astropy.io.fits import getdata
from pathlib import Path
import healpy as hp
import sqlalchemy, json, os, sys, glob
from qa_ga_sim import (
    export_results
)

confg_qa = "qa_ga_sim.json"

with open(confg_qa) as fstream:
    param_qa = json.load(fstream)

globals().update(param_qa)

with open(pars_simulations_file) as fstream:
    param = json.load(fstream)

os.system('jupyter nbconvert --execute --to html --EmbedImagesPreprocessor.embed_images=True qa_ga_sim.ipynb')

# export_results(param['export_path'], param['results_path'], param['copy_html_path'])

os.makedirs(os.path.dirname(export_file), exist_ok=True)

os.system('cp qa_ga_sim.html ' + export_file)