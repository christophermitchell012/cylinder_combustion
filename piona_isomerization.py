"""Model a restricted C5 isomerization pass and design a paired GC-VUV challenge.

Run: python -m pip install cantera==3.2.0
     python piona_isomerization.py --help

Fractions describe mass throughout. P, I, O, N and A denote normal paraffins,
isoparaffins, olefins, naphthenes and aromatics. Nothing is normalized to the
sum of the reported classes; oxygenates and other matrix material may exist.
"""
import argparse
import csv
import math
from pathlib import Path

import cantera as ct

CLASSES = ('P', 'I', 'O', 'N', 'A')
NAMES = ('C5H12,n-pentane', 'C5H12,i-pentane')


def bounded(value, low, high, label):
    if not math.isfinite(value) or not low <= value <= high:
        raise ValueError(f'{label} must be finite and within [{low}, {high}]')
    return value


def equilibrium_curve(temperatures):
    """Evaluate n-pentane ⇌ 2-methylbutane in an ideal vapor at fixed T."""
    species = {s.name: s for s in ct.Species.list_from_file('nasa_gas.yaml')}
    gas = ct.Solution(thermo='ideal-gas', species=[species[n] for n in NAMES])
    for kelvin in temperatures:
        bounded(kelvin, gas.min_temp, gas.max_temp, 'temperature K')
        gas.TP = kelvin, ct.one_atm
        # At equilibrium, mu_n = mu_i. Equal stoichiometric mole counts cancel
        # pressure, so ln(x_i/x_n) = (G_n° - G_i°)/(RT). Equal MW makes x_i=w_i.
        g_normal, g_iso = gas.standard_gibbs_RT
        log_ratio = g_normal - g_iso
        iso = 1 / (1 + math.exp(-log_ratio))
        yield kelvin, iso


def outlet(feed_iso, eq_iso, approach):
    # An independently specified fraction of the distance to equilibrium.
    # This scenario parameter is neither measured conversion nor a kinetic law.
    return feed_iso + approach * (eq_iso - feed_iso)


def paired_delta(iso_a, iso_b, binary_fraction):
    # Equal fractions of one unchanged matrix cancel even if its assay is unknown.
    shift = 100 * binary_fraction * (iso_b - iso_a)
    return dict(zip(CLASSES, (-shift, shift, 0.0, 0.0, 0.0)))


def compare(path):
    """Compare observed B-A with independent, as-weighed binary standards."""
    required = {'sample', 'normal_g', 'iso_g', 'matrix_g', *CLASSES}
    samples = {}
    with path.open(newline='', encoding='utf-8-sig') as stream:
        reader = csv.DictReader(stream)
        if not reader.fieldnames or set(reader.fieldnames) != required or len(reader.fieldnames) != len(required):
            raise ValueError('CSV header must be sample,normal_g,iso_g,matrix_g,P,I,O,N,A')
        for row in reader:
            label = row['sample']
            if label not in ('A', 'B') or label in samples or None in row:
                raise ValueError('expected exactly one A and one B, with one value per column')
            def number(key, upper):
                return bounded(float(row[key]), 0, upper, f'{label}/{key}')
            grams = [number(k, 1e100) for k in ('normal_g', 'iso_g', 'matrix_g')]
            total = sum(grams)
            if total <= 0:
                raise ValueError('empty sample')
            assay = {k: number(k, 100) for k in CLASSES}
            if sum(assay.values()) > 100.05:
                raise ValueError('class sum exceeds 100 mass % beyond rounding')
            samples[label] = ([g/total for g in grams], assay)
    if set(samples) != {'A', 'B'}:
        raise ValueError('both samples A and B are required')
    (a, observed_a), (b, observed_b) = samples['A'], samples['B']
    # Matrix mismatch bounds the unknown class correction: its fraction in
    # the common matrix lies between zero and one. Keep that uncertainty visible.
    matrix_term = 100*(b[2]-a[2])
    print('\nAs-weighed comparison, percentage points (B minus A)')
    print('class  expected interval    observed   residual to nearest bound')
    for index, group in enumerate(CLASSES):
        spike = 100*(b[index]-a[index]) if index < 2 else 0.0
        low, high = sorted((spike, spike+matrix_term))
        obs = observed_b[group]-observed_a[group]
        residual = obs - min(max(obs, low), high)
        print(f'{group:>3}   [{low:8.3f},{high:8.3f}]  {obs:9.3f}  {residual:+9.3f}')


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--temperature-k', type=float, nargs='+', default=[423.15, 453.15, 483.15])
    parser.add_argument('--feed-iso', type=float, default=0.20, help='feed iso mass fraction within binary C5')
    parser.add_argument('--approach', type=float, default=0.80, help='fraction of distance to equilibrium')
    parser.add_argument('--binary-fraction', type=float, default=0.10, help='binary mass fraction in final blend')
    parser.add_argument('--compare', type=Path, metavar='CSV', help='as-weighed A/B masses and reported class mass percentages')
    args = parser.parse_args()
    try:
        for name in ('feed_iso', 'approach', 'binary_fraction'):
            bounded(getattr(args, name), 0, 1, name)
        rows = list(equilibrium_curve(args.temperature_k))
        print('     K    eq I/C5   outlet I/C5   delta P pp   delta I pp')
        for temp, eq in rows:
            product = outlet(args.feed_iso, eq, args.approach)
            delta = paired_delta(args.feed_iso, product, args.binary_fraction)
            print(f'{temp:7.2f}  {eq:9.4f}  {product:12.4f}  {delta["P"]:11.3f}  {delta["I"]:11.3f}')
        if args.compare:
            compare(args.compare)
    except (ValueError, TypeError, OSError, csv.Error, ct.CanteraError) as exc:
        parser.exit(2, f'Input or model error: {exc}\n')


if __name__ == '__main__':
    main()
