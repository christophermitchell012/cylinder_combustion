"""Radial temperature and pressure snapshots with a grid-refinement comparison.

Run with Python 3. Install numpy, scipy, matplotlib, pillow, and cantera.
The generated image is saved in the current working directory.
"""

import numpy as np
import cantera as ct
import matplotlib.pyplot as plt
from scipy.special import expit
from matplotlib.ticker import MaxNLocator
# Use GRI species for equilibrium products and NASA thermodynamic data for iso-octane.
sp=ct.Species.list_from_file('gri30.yaml')
fuel=next(s for s in ct.Species.list_from_file('nasa_gas.yaml') if s.name=='C8H18,isooctane')
gas=ct.Solution(thermo='ideal-gas',species=sp+[fuel])
# Compressed fresh charge at 700 K and 15 bar. Cantera expects pressure in Pa.
gas.TP=700,1.5e6
gas.set_equivalence_ratio(1,{fuel.name:1},{'O2':1,'N2':3.76})
rho,a,g=gas.density,gas.sound_speed,gas.cp_mass/gas.cv_mass
# Hold enthalpy and pressure fixed to estimate the fully burned temperature and expansion.
gas.equilibrate('HP')
Tb,ratio=gas.T,gas.density/rho
R=.043
def run(n):
    dr=R/n
    # Pressure lives at annular cell centers and radial velocity lives at faces.
    r=(np.arange(n)+.5)*dr
    face=np.arange(n+1)*dr
    # Subdivide each output interval to keep the acoustic Courant number at or below 0.4.
    dt=25e-6/np.ceil(25e-6*a/(.4*dr))
    p=np.zeros(n)
    u=np.zeros(n+1)
    # Wiebe burned-mass fraction reaches 99 percent 50 crank degrees after ignition.
    def burned(t):
        # Prescribe a 6 µs ignition ramp and a smooth front expanding at 18 m/s.
        return (1-np.exp(-(t/6e-6)**2))*expit((.0006+18*t-r)/.00012)
    rows=[]
    for i in range(round(125e-6/dt)+1):
        t=i*dt
        if i%round(25e-6/dt)==0:
            # Annular volumes scale with radius at constant cell width and chamber depth.
            avg=np.average(p,weights=r)
            rows.append((t,700+(Tb-700)*burned(t),p.copy(),avg))
        # Advance momentum at interior faces. Zero velocity at both endpoints enforces symmetry and a rigid wall.
        u[1:-1]-=dt*np.diff(p)/(rho*dr)
        # Cylindrical divergence includes the larger area of the outer cell face.
        p-=rho*a*a*dt*np.diff(face*u)/(r*dr)
        # Newly burned volume acts as an acoustic expansion source through the equilibrium density ratio.
        p+=rho*a*a*(1-ratio)*(burned(t+dt)-burned(t))
    return r,rows
r,rows=run(600)
# Repeat with half-width cells to estimate sensitivity to radial resolution.
rf,fine=run(1200)
print(a,Tb,ratio,[(round(t*1e6),round(m),round(max(abs(p-m)))) for t,T,p,m in rows])
print('grid diff Pa',max(np.max(abs(np.interp(r,rf,q[2])-z[2])) for z,q in zip(rows,fine)))
fig,axes=plt.subplots(2,1,figsize=(9,7),layout='constrained')
for ax in axes: ax.grid(alpha=.22);ax.set_xlabel('Distance from plug (mm)')
for color,(t,T,p,m) in zip(plt.cm.viridis(np.linspace(.09,.87,6)),rows):
    label=f'{t*1e6:.0f} µs'
    axes[0].plot(r*1000,T,color=color,label=label,lw=2)
    # Subtract the instantaneous chamber mean to isolate the traveling pressure disturbance.
    axes[1].plot(r*1000,(p-m)/1000,color=color,lw=1.7)
axes[0].set(xlim=(0,4),ylabel='Temperature (K)',title='Flame travel • first 4 mm from the plug')
axes[0].legend(ncol=3)
axes[1].set(xlim=(0,43),ylabel='Pressure minus chamber mean (kPa)',title='Pressure wave • full 43 mm radius, reflecting wall at right')
fig.suptitle('Center-ignited cylinder  |  snapshots every 25 µs',fontsize=16)
fig.savefig('cylinder-propagation.png',dpi=175)
