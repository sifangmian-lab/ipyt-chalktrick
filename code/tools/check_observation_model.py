"""Physical invariants for the finite-footprint observation model."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from observation_model import footprint_kernel,profile_to_pixel_mass,render_mass,camera_response

def main():
    x=np.linspace(.5,99.5,100)*1e-6
    i=2+np.sin(np.arange(100)*.1)
    grid,m=profile_to_pixel_mass(x,i,9.3)
    original=np.sum(i)*(x[1]-x[0])
    assert abs(m.sum()/original-1)<1e-12
    k=footprint_kernel(27,1.5,.31)
    _,_,d=render_mass(m,k,grid[0],.9)
    assert abs(d.sum()/original-1)<1e-11
    # Linear superposition before camera response, including PSF.
    _,_,d1=render_mass(m*.3,k,grid[0],.9)
    _,_,d2=render_mass(m*.7,k,grid[0],.9)
    assert np.allclose(d1+d2,d,rtol=1e-10,atol=1e-16)
    # Unit point footprint and zero blur recover the original source exactly.
    _,_,point=render_mass(m,np.ones((1,1)),grid[0],0)
    assert np.allclose(point[0],m,atol=1e-16)
    # A finite tilted footprint adds x variance despite identical force/source.
    impulse=np.ones(1)
    x0,_,upright=render_mass(impulse,footprint_kernel(27,1.5,0),0,0)
    xt,_,tilted=render_mass(impulse,k,0,0)
    variance=lambda xx,dd:np.sum(xx*xx*dd.sum(axis=0))/dd.sum()
    assert variance(xt,tilted)>variance(x0,upright)+10
    # Nonlinear camera response is applied per pixel before transverse integration.
    a=np.array([[1.,0.],[1.,0.]])
    nonlinear=camera_response(a,gain=2,mode='saturating').sum(axis=0)
    wrong=camera_response(a.sum(axis=0),gain=2,mode='saturating')
    assert not np.allclose(nonlinear,wrong)
    assert np.all(camera_response(a,gain=2,mode='clipped')<=1)
    print('PASS: source mass, footprint/PSF mass, superposition, point limit, geometric broadening, camera order.')

if __name__=='__main__':main()
