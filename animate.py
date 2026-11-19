"""Top-down animation of ignition and acoustic propagation in a fixed disk.

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
    dt=2.5e-6/np.ceil(2.5e-6*a/(.4*dr))
    p=np.zeros(n)
    u=np.zeros(n+1)
    # Wiebe burned-mass fraction reaches 99 percent 50 crank degrees after ignition.
    def burned(t):
        # Prescribe a 6 µs ignition ramp and a smooth front expanding at 18 m/s.
        return (1-np.exp(-(t/6e-6)**2))*expit((.0006+18*t-r)/.00012)
    rows=[]
    for i in range(round(125e-6/dt)+1):
        t=i*dt
        if i%round(2.5e-6/dt)==0:
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
from PIL import Image
r,rows=run(600)
fig,axes=plt.subplots(1,3,figsize=(10,4),layout='constrained')
frames=[]
for t,T,p,m in rows:
    for ax,L,field,cmap,low,high,title in zip(axes,[43,4,43],[T,T,(p-m)/1000],['inferno','inferno','RdBu_r'],[700,700,-25],[2551,2551,25],['Temperature • full bore','Temperature • plug close-up','Pressure • full bore']):
        if t==0:
            x=np.linspace(-L,L,251)
            # Map each image pixel to its radial distance from the center-mounted plug.
            ax.rad=np.hypot(x[:,None],x)
            ax.im=ax.imshow(np.zeros_like(ax.rad),extent=(-L,L,-L,L),origin='lower',cmap=cmap,vmin=low,vmax=high)
            ax.set_title(title,fontsize=11)
            ax.set_xlabel('mm')
            if L==43: ax.add_patch(plt.Circle((0,0),43,fill=False,color='gray',lw=.8))
            fig.colorbar(ax.im,ax=ax,orientation='horizontal',shrink=.86,label='Temperature (K)' if low==700 else 'Local minus mean (kPa)',extend='both' if low<0 else 'neither')
        # Interpolate the radial solution in meters and mask pixels outside the bore.
        ax.im.set_data(np.ma.masked_where(ax.rad>43,np.interp(ax.rad/1000,r,field)))
    fig.suptitle(f'Center ignition    |    {t*1e6:05.1f} µs    |    Mean pressure {15+m/1e5:.3f} bar')
    fig.canvas.draw()
    frames.append(Image.fromarray(np.asarray(fig.canvas.buffer_rgba())[:,:,:3].copy()))
path='cylinder-animation.gif'
# Hold the final frame briefly before replaying the simulation.
frames[0].save(path,save_all=True,append_images=frames[1:],duration=[100]*50+[900],loop=0)
im=Image.open(path)
print(im.n_frames,im.size)
