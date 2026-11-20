"""Compare the compiled CHEMKIN-II driver with the existing Cantera model."""
import csv
import io
from pathlib import Path
import subprocess
import sys

import numpy as np

root = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(root))
from piona_isomerization import equilibrium_curve, outlet, paired_delta


def compare(executable, temperatures):
    run = subprocess.run([executable], input=''.join(f'{t:.12f}\n' for t in temperatures),
                         text=True, capture_output=True, check=True)
    rows = list(csv.DictReader(io.StringIO(run.stdout)))
    if len(rows) != len(temperatures):
        raise RuntimeError(f'Incomplete CHEMKIN run: {run.stderr}')
    differences = []
    for row, (temp, eqIso) in zip(rows, equilibrium_curve(temperatures)):
        if abs(float(row['temperature_K'])-temp) > 1e-10:
            raise ValueError('Temperature mismatch')
        product = outlet(0.2, eqIso, 0.8)
        delta = paired_delta(0.2, product, 0.1)
        reference = [eqIso, product, delta['P'], delta['I']]
        result = [float(row[key]) for key in ('equilibrium_iso','outlet_iso','delta_P_pp','delta_I_pp')]
        differences.append(np.abs(np.subtract(reference,result)))
    peak = np.max(differences, axis=0)
    print(f'{len(rows)} temperatures, {temperatures[0]:g} to {temperatures[-1]:g} K')
    print('Maximum absolute differences (eq fraction, outlet fraction, P pp, I pp):', peak)
    if not np.all(peak < 1e-9):
        raise RuntimeError('Implementations differ by more than the comparison tolerance')


if __name__ == '__main__':
    executable = Path(sys.argv[1] if len(sys.argv)>1 else './piona').resolve()
    compare(executable, np.linspace(423.15,483.15,601))
    # Test both NASA7 polynomial intervals and temperatures next to the join.
    compare(executable, np.unique(np.r_[298.15, np.arange(299,5001),999.999999,1000.000001]))
