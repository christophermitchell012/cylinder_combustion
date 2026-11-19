"""Vertical section through a moving piston with prescribed two-zone combustion.

Run with Python 3. Install numpy, scipy, matplotlib, pillow, and cantera.
The generated image is saved in the current working directory.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from PIL import Image
# Square 86 mm bore/stroke, 143 mm rod, and 10:1 compression ratio. Geometry uses meters.
D=.086; A=np.pi*D**2/4; r=D/2; rod=.143; clearance=D/9; gamma=1.32
def height(d):
 # Slider-crank displacement is measured from top dead center.
 t=np.deg2rad(d)
 return clearance+r*(1-np.cos(t))+rod-np.sqrt(rod**2-(r*np.sin(t))**2)
def volume(d):return A*height(d)
# Wiebe burned-mass fraction reaches 99 percent 50 crank degrees after ignition.
def burn(d):return 1-np.exp(-4.605*((d+15)/50)**3)
def rhs(d,y):
 # Differentiate chamber volume with respect to crank angle in degrees.
 V=volume(d); dv=(volume(d+.0001)-volume(d-.0001))/.0002
 # Analytic derivative of the prescribed burned-mass fraction per crank degree.
 dx=np.exp(-4.605*((d+15)/50)**3)*3*4.605*(d+15)**2/50**3
 # Closed-cylinder first law balances heat release against piston work at constant gamma.
 return ((gamma-1)*mass*2.6e6*dx-gamma*y[0]*dv)/V
# Estimate trapped mass from the ideal gas law at the ignition angle.
p0=1.5e6; T0=700; mass=p0*volume(-15)/(287*T0)
angles=np.arange(-15,91)
ps=solve_ivp(rhs,(-15,90),[p0],rtol=1e-10,atol=.001,t_eval=angles).y[0]
fig,(ax,trace)=plt.subplots(1,2,figsize=(9,4.6),gridspec_kw={'width_ratios':[1.18,1]},layout='constrained')
X,Z=np.meshgrid(np.linspace(-43,43,260),np.linspace(0,65,200))
frames=[]
for d,p in zip(angles,ps):
 h=height(d)*1000; V=volume(d); x=burn(d); Tu=T0*(p/p0)**((gamma-1)/gamma)
 # Isentropic compression of the fresh gas sets its volume at the shared cylinder pressure.
 unburned=(1-x)*mass*287*Tu/p; f=np.clip(1-unburned/V,0,1)
 # Recover burned-zone temperature from its mass, volume, and the common gas constant.
 Tb=p*V*f/(max(x,1e-12)*mass*287) if x>1e-10 else Tu
 def fraction(s):
  # Clip the flame hemisphere at the piston depth and at the 43 mm bore radius.
  b=min(s,h); k=min(np.sqrt(max(s*s-43**2,0)),b)
  return (43**2*k+s*s*(b-k)-(b**3-k**3)/3)/(43**2*h)
 # Find the hemisphere radius enclosing the required burned volume in the current chamber.
 s=brentq(lambda s:fraction(s)-f,0,np.hypot(43,h)+.01) if 0<f<1-1e-12 else (0 if f==0 else np.hypot(43,h))
 ax.clear();trace.clear()
 # Classify each pixel in the vertical slice by its distance from the spark.
 section=np.where(np.hypot(X,Z)<=s,Tb,Tu)
 im=ax.pcolormesh(X,Z,np.ma.masked_where(Z>h,section),cmap='inferno',vmin=650,vmax=3250,shading='auto')
 # Draw the piston crown at the slider-crank height for this frame.
 ax.fill_between([-43,43],h,h+7,color='#566674');ax.plot([-43,43],[h,h],color='black',lw=2)
 ax.plot([-43,-43,43,43],[63,0,0,63],color='#334155',lw=3)
 ax.plot(0,0,'v',color='#13cce0',ms=7);ax.set(xlim=(-47,47),ylim=(65,-4),aspect='equal',xlabel='Distance from bore axis (mm)',ylabel='Depth below head (mm)',title='Gas temperature • centerline section')
 if d==-15:fig.colorbar(im,ax=ax,orientation='horizontal',shrink=.83,label='Gas temperature (K)')
 trace.plot(angles[angles<=d],ps[angles<=d]/1e5,color='#176789',lw=2.5);trace.plot(d,p/1e5,'o',color='#176789');trace.axvline(0,color='gray',ls=':')
 trace.set(xlim=(-15,90),ylim=(0,np.ceil(max(ps)/1e6)*10),xlabel='Crank angle (degrees from TDC)',ylabel='Absolute pressure (bar)',title='Chamber-average pressure');trace.grid(alpha=.2)
 fig.suptitle(f'{d:+d}° from TDC     |     {p/1e5:.1f} bar     |     {100*x:.0f}% mass burned',fontsize=14)
 fig.canvas.draw();frames.append(Image.fromarray(np.asarray(fig.canvas.buffer_rgba())[:,:,:3].copy()))
out='cylinder-cross-section.gif'
# Hold the final frame briefly before replaying the simulation.
frames[0].save(out,save_all=True,append_images=frames[1:],duration=[70]*105+[800],loop=0)
print(Image.open(out).n_frames,ps.max()/1e5, angles[np.argmax(ps)],height(90)*1000)
