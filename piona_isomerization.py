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

classCodes = ('P', 'I', 'O', 'N', 'A')
pentaneNames = ('C5H12,n-pentane', 'C5H12,i-pentane')


def bounded(value, low, high, label):
    if not math.isfinite(value) or not low <= value <= high:
        raise ValueError(f'{label} must be finite and within [{low}, {high}]')
    return value


def equilibrium_curve(temperatures):
    """Evaluate n-pentane ⇌ 2-methylbutane in an ideal vapor at fixed T."""
    species = {s.name: s for s in ct.Species.list_from_file('nasa_gas.yaml')}
    gas = ct.Solution(thermo='ideal-gas', species=[species[n] for n in pentaneNames])
    for kelvin in temperatures:
        bounded(kelvin, gas.min_temp, gas.max_temp, 'temperature K')
        gas.TP = kelvin, ct.one_atm
        # One pentane molecule becomes one isopentane molecule, so pressure drops
        # out of the ratio. Both weigh the same, so mole and mass fractions agree.
        normalGibbsRT, isoGibbsRT = gas.standard_gibbs_RT
        logIsoRatio = normalGibbsRT - isoGibbsRT
        iso = 1 / (1 + math.exp(-logIsoRatio))
        yield kelvin, iso


def outlet(feedIso, eqIso, approach):
    # At 0 the feed passes through unchanged. At 1 it reaches equilibrium.
    # Pick a value in between to model a bed that only gets partway there.
    return feedIso + approach * (eqIso - feedIso)


def paired_delta(isoA, isoB, binaryFraction):
    # Both blends get the same share of base gasoline. Subtracting A from B
    # cancels that gasoline out, even without a seperate assay.
    shift = 100 * binaryFraction * (isoB - isoA)
    return dict(zip(classCodes, (-shift, shift, 0.0, 0.0, 0.0)))


def compare(path):
    """Compare observed B-A with independent, as-weighed binary standards."""
    required = {'sample', 'normal_g', 'iso_g', 'matrix_g', *classCodes}
    samples = {}
    with path.open(newline='', encoding='utf-8-sig') as stream:
        reader = csv.DictReader(stream)
        if not reader.fieldnames or set(reader.fieldnames) != required or len(reader.fieldnames) != len(required):
            raise ValueError('CSV header must be sample,normal_g,iso_g,matrix_g,P,I,O,N,A')
        for row in reader:
            label = row['sample']
            if label not in ('A', 'B') or label in samples or None in row or None in row.values():
                raise ValueError('expected exactly one A and one B, with one value per column')
            def number(key, upper):
                return bounded(float(row[key]), 0, upper, f'{label}/{key}')
            grams = [number(k, 1e100) for k in ('normal_g', 'iso_g', 'matrix_g')]
            total = sum(grams)
            if total <= 0:
                raise ValueError('empty sample')
            assay = {k: number(k, 100) for k in classCodes}
            if sum(assay.values()) > 100.05:
                raise ValueError('class sum exceeds 100 mass % beyond rounding')
            samples[label] = ([g/total for g in grams], assay)
    if set(samples) != {'A', 'B'}:
        raise ValueError('both samples A and B are required')
    (weightsA, assayA), (weightsB, assayB) = samples['A'], samples['B']
    # If the matrix shares differ, some of the shift could come from the gasoline.
    # Bound that contribution by allowing each class to make up 0–100% of it.
    matrixShift = 100 * (weightsB[2] - weightsA[2])
    print('\nAs-weighed comparison, percentage points (B minus A)')
    print('class  expected interval    observed   residual to nearest bound')
    for index, group in enumerate(classCodes):
        spike = 100 * (weightsB[index] - weightsA[index]) if index < 2 else 0.0
        low, high = sorted((spike, spike + matrixShift))
        observedShift = assayB[group] - assayA[group]
        residual = observedShift - min(max(observedShift, low), high)
        print(f'{group:>3}   [{low:8.3f},{high:8.3f}]  {observedShift:9.3f}  {residual:+9.3f}')


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--temperature-k', dest='tempK', type=float, nargs='+', default=[423.15, 453.15, 483.15])
    parser.add_argument('--feed-iso', dest='feedIso', type=float, default=0.20, help='feed iso mass fraction within binary C5')
    parser.add_argument('--approach', type=float, default=0.80, help='fraction of distance to equilibrium')
    parser.add_argument('--binary-fraction', dest='binaryFraction', type=float, default=0.10, help='binary mass fraction in final blend')
    parser.add_argument('--compare', type=Path, metavar='CSV', help='as-weighed A/B masses and reported class mass percentages')
    args = parser.parse_args()
    try:
        for name in ('feedIso', 'approach', 'binaryFraction'):
            bounded(getattr(args, name), 0, 1, name)
        rows = list(equilibrium_curve(args.tempK))
        print('     K    eq I/C5   outlet I/C5   delta P pp   delta I pp')
        for temp, eqIso in rows:
            outletIso = outlet(args.feedIso, eqIso, args.approach)
            delta = paired_delta(args.feedIso, outletIso, args.binaryFraction)
            print(f'{temp:7.2f}  {eqIso:9.4f}  {outletIso:12.4f}  {delta["P"]:11.3f}  {delta["I"]:11.3f}')
        if args.compare:
            compare(args.compare)
    except (ValueError, TypeError, OSError, csv.Error, ct.CanteraError) as exc:
        parser.exit(2, f'Input or model error: {exc}\n')


if __name__ == '__main__':
    main()
