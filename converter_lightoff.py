"""Compare cold-start CO using a CatMAP Pt(111) rate table and a lumped brick."""
import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import cumulative_trapezoid, solve_ivp
from scipy.interpolate import RegularGridInterpolator

DATA = Path(__file__).resolve().parent / 'converter'
AMBIENT_K = 293.15
BRICK_J_K = 520.0           # Effective thermally active brick heat capacity
SITES_MOL = 0.001          # Accessible sites, an illustrative inventory
MW_CO = 28.01              # g/mol
MW_EXHAUST = 0.029         # kg/mol


def charge_mass():
    """Reproduce stroke.py's ignition-state charge without running its animation."""
    radius, rod = 0.043, 0.143
    angle = np.deg2rad(-15)
    height = 0.086/9 + radius*(1-np.cos(angle)) + rod - np.sqrt(
        rod**2 - (radius*np.sin(angle))**2)
    return 1.5e6 * np.pi * radius**2 * height / (287*700)


def inlet(t):
    """Synthetic fast idle, four cylinders, with time in seconds after ignition."""
    rpm = 1050 + 450*np.exp(-t/30)
    mass_flow = charge_mass()*4*rpm/120
    co_fraction = 0.002 + 0.012*np.exp(-t/18)
    gas_k = 560 + 320*(1-np.exp(-t/15))
    return mass_flow, co_fraction, gas_k


def brick_temperature(time, initial_k, heater_w=0, heater_end=0):
    """Uniform brick, effective gas-to-solid transfer, and ambient heat loss."""
    def rhs(t, temp):
        flow, _, gas_k = inlet(t)
        heat = 0.85*flow*1100*(gas_k-temp) - 2*(temp-AMBIENT_K)
        return (heat + (heater_w if t < heater_end else 0))/BRICK_J_K

    return solve_ivp(rhs, (0, time[-1]), [initial_k], t_eval=time,
                     max_step=0.5, rtol=1e-8, atol=1e-7).y[0]


def rate_lookup():
    table = np.loadtxt(DATA/'rates.csv', delimiter=',', skiprows=1)
    temp, pressure = (np.unique(table[:, i]) for i in (0, 1))
    rates = table[:, 2].reshape(len(temp), len(pressure))
    # Interpolate log turnover because the rate spans orders of magnitude.
    return RegularGridInterpolator((temp, pressure), np.log(np.maximum(rates, 1e-100)))


def simulate(time, initial_k, lookup, heater_w=0, heater_end=0, segments=80):
    temp = brick_temperature(time, initial_k, heater_w, heater_end)
    mass_flow, x_co, _ = inlet(time)
    total_in = mass_flow / MW_EXHAUST
    co_in = total_in*x_co
    co = co_in.copy()
    for _ in range(segments):
        # The inlet O2:CO ratio is 1:2. Reaction preserves this ratio and
        # contracts total molar flow by half a mole for every mole CO oxidized.
        total = total_in - 0.5*(co_in-co)
        reactive_bar = 1.5*co/total  # Total exhaust pressure is 1 bar.
        active = (reactive_bar >= 1e-7) & (temp >= 300)
        points = np.column_stack((temp[active], np.log10(reactive_bar[active])))
        capacity = np.exp(lookup(points))*SITES_MOL/segments
        co[active] -= np.minimum(co[active], capacity)
    flow_g_s = co*MW_CO
    cumulative = cumulative_trapezoid(flow_g_s, time, initial=0)
    return temp, flow_g_s, cumulative


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=Path('.'))
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    lookup = rate_lookup()
    time = np.arange(0, 60.01, 0.25)
    mass, fraction, _ = inlet(time)
    inlet_g_s = mass/MW_EXHAUST*fraction*MW_CO
    fig, axes = plt.subplots(2, sharex=True, constrained_layout=True)
    axes[0].plot(time, inlet_g_s*1000, ':k', label='Engine out')
    columns = [time, inlet_g_s]
    for name, offset in [('Cold brick', 0), ('Preheated +100 °C', 100)]:
        temp, flow, grams = simulate(time, AMBIENT_K+offset, lookup)
        axes[0].plot(time, flow*1000, label=name)
        axes[1].plot(time, grams, label=name)
        columns.extend([temp, flow, grams])
        print(f'{name}: {grams[-1]:.3f} g tailpipe CO in 60 s')
    axes[0].set_title('Gasoline cold start | illustrative four-cylinder fast idle')
    axes[0].set_ylabel('CO flow (mg/s)')
    axes[1].set_ylabel('Tailpipe CO (g)')
    axes[1].set_xlabel('Time after ignition (s)')
    axes[0].legend()
    fig.savefig(args.output_dir/'converter-lightoff.png', dpi=160)
    plt.close(fig)
    np.savetxt(args.output_dir/'converter-lightoff.csv', np.column_stack(columns),
               delimiter=',', comments='', header='time_s,inlet_g_s,cold_K,cold_g_s,cold_g,preheat_K,preheat_g_s,preheat_g')
    # Only Base and Case 1 have no secondary-air injection. Match their
    # heater schedule and observation window, retaining our synthetic inlet.
    with (DATA/'bench.csv').open() as stream:
        cases = list(csv.DictReader(stream))
    for case in cases[:2]:
        watts = 1000*float(case['heater_kw'])
        pre_seconds = -float(case['heater_start_s'])
        # Before ignition there is no exhaust. Solve the heat-loss ODE exactly.
        offset = watts/2*(-np.expm1(-2*pre_seconds/BRICK_J_K))
        clock = np.arange(0, float(case['window_s'])+0.25, 0.25)
        _, _, grams = simulate(clock, AMBIENT_K+offset, lookup, watts,
                               float(case['heater_end_s']))
        observed = float(case['measured_tailpipe_co_g'])
        print(f"{case['case']} at 500 s: model {grams[-1]:.2f} g, "
              f"bench {observed:.2f} g, difference {grams[-1]-observed:+.2f} g")


if __name__ == '__main__':
    main()
