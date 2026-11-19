"""Rebuild rates.csv with CatMAP in a clean directory to avoid stale solver caches."""
from pathlib import Path
import os
import shutil
import tempfile

from catmap import ReactionModel
import numpy as np

root = Path(__file__).resolve().parent
with tempfile.TemporaryDirectory() as work:
    for name in ('pt.mkm', 'energies.txt'):
        shutil.copy(root / name, work)
    original = Path.cwd()
    try:
        os.chdir(work)
        model = ReactionModel(setup_file='pt.mkm')
        model.output_variables += ['production_rate']
        model.run()
        index = model.output_labels['production_rate'].index('CO2_g')
        points = {tuple(point): float(rate[index]) for point, rate in model.production_rate_map}
        rows = []
        for temp in np.linspace(*model.descriptor_ranges[0], model.resolution[0]):
            for pressure in np.linspace(*model.descriptor_ranges[1], model.resolution[1]):
                rate = points[(temp, pressure)]
                if not np.isfinite(rate) or rate < 0:
                    raise ValueError(f'Invalid turnover at {temp} K, log pressure {pressure}')
                rows.append((temp, pressure, rate))
        np.savetxt('rates.csv', rows, fmt='%.12g', delimiter=',',
                   header='temperature_K,log10_reactive_pressure_bar,co2_per_site_s', comments='')
        shutil.copy('rates.csv', root / 'rates.csv')
    finally:
        os.chdir(original)
