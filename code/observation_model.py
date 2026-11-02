"""Finite-footprint forward observation model, in image-pixel coordinates.

Input masses are per-bin integrated deposit proxies, NOT grayscale or force.
K is a nonnegative unit-sum pixel mass kernel. Convolution and optical PSF
preserve total proxy mass on padded grids; camera response generally does not.
Footprint, length scale, optical gain and blur require independent calibration.
"""
from __future__ import annotations
import numpy as np
from scipy import ndimage, signal


def footprint_kernel(half_height_px=27, sigma_x_px=1.5, tilt_dx_dy=.31):
    """Stationary sheared ridge with semicircular y envelope, sum(K)=1.

    tilt is dx/dy, not a chalk angle. This simple footprint has no grain model.
    It is a phenomenological observation kernel, not a measured contact patch.
    """
    if half_height_px < 0 or sigma_x_px < 0 or not np.isfinite(tilt_dx_dy):
        raise ValueError('Invalid footprint geometry')
    if half_height_px == 0 and sigma_x_px == 0:
        return np.ones((1, 1))
    if sigma_x_px <= 0:
        raise ValueError('A finite footprint needs positive sigma_x_px')
    h=int(np.ceil(half_height_px))
    w=int(np.ceil(abs(tilt_dx_dy)*half_height_px+5*sigma_x_px))
    y=np.arange(-h,h+1)[:,None];x=np.arange(-w,w+1)[None,:]
    envelope=np.sqrt(np.maximum(1-(y/max(half_height_px,1e-12))**2,0)) if h else np.ones((1,1))
    kernel=envelope*np.exp(-.5*((x-tilt_dx_dy*y)/sigma_x_px)**2)
    if kernel.sum() <= 0:
        raise ValueError('Empty footprint')
    return kernel/kernel.sum()


def profile_to_pixel_mass(x, intensity, span_px):
    """Uniform-bin I(x) -> pixel masses by linear assignment, preserving sum(I dx).

    span_px is an explicitly chosen image scale. It is not inferred from data.
    Model x can be metres; output image coordinates are pixels centred at zero.
    """
    x=np.asarray(x,dtype=float);i=np.asarray(intensity,dtype=float)
    if x.ndim!=1 or len(x)<2 or x.shape!=i.shape or span_px<=0:
        raise ValueError('Need >=2 bins and positive span')
    dx=np.diff(x)
    if not np.all(dx>0) or not np.allclose(dx,dx[0],rtol=2e-4):
        raise ValueError('Input bins must be uniformly increasing')
    if not np.all(np.isfinite(i)) or np.any(i<0):
        raise ValueError('Deposited proxy mass must be finite and nonnegative')
    source_mass=i*dx[0]
    physical_span=(x[-1]-x[0])+dx[0]
    coordinates=((x-(x[0]+x[-1])/2)/physical_span)*span_px
    radius=int(np.ceil(span_px/2))+2
    grid=np.arange(-radius,radius+1,dtype=float)
    mass=np.zeros(len(grid))
    positions=coordinates+radius
    left=np.floor(positions).astype(int);fraction=positions-left
    np.add.at(mass,left,source_mass*(1-fraction))
    np.add.at(mass,left+1,source_mass*fraction)
    return grid,mass


def render_mass(source_mass, kernel, source_x0=0., psf_sigma_px=.7):
    """Full convolution followed by padded Gaussian PSF. Returns x,y,mass2d."""
    source=np.asarray(source_mass,dtype=float);k=np.asarray(kernel,dtype=float)
    if source.ndim!=1 or k.ndim!=2 or np.any(source<0) or np.any(k<0):
        raise ValueError('Nonnegative 1D source and 2D kernel required')
    if not np.isclose(k.sum(),1,rtol=1e-10,atol=1e-12) or psf_sigma_px<0:
        raise ValueError('Kernel must have unit mass and PSF must be nonnegative')
    field=signal.fftconvolve(source[None,:],k,mode='full')
    # Roundoff can generate tiny negative values after FFT convolution.
    field=np.maximum(field,0)
    pad=int(np.ceil(4*psf_sigma_px))+1 if psf_sigma_px else 0
    if pad:
        field=np.pad(field,pad)
        field=ndimage.gaussian_filter(field,psf_sigma_px,mode='constant',truncate=4)
    x=np.arange(field.shape[1])+source_x0-(k.shape[1]-1)/2-pad
    y=np.arange(field.shape[0])-(k.shape[0]-1)/2-pad
    return x,y,field


def camera_response(mass2d, gain=1., background=0., mode='linear', white=1.):
    """Illustrative camera map; nonlinear modes are uncalibrated scenarios."""
    z=np.asarray(mass2d,dtype=float)*gain
    if gain<0 or white<=background:
        raise ValueError('Invalid camera gain/range')
    if mode=='linear':
        return background+z
    if mode=='saturating':
        return background+(white-background)*(-np.expm1(-z/(white-background)))
    if mode=='clipped':
        return np.minimum(background+z,white)
    raise ValueError('Unknown camera mode')


def camera_image_on_crop(x,y,field,shape=(160,96),x_shift=0.,y_center=80.,gain=1.):
    """Sample a forward field on the unrotated crop grid, preserving original x."""
    rows,cols=np.indices(shape,dtype=float)
    coordinates=np.array([rows-y_center-y[0],cols-shape[1]//2-x_shift-x[0]])
    return gain*ndimage.map_coordinates(field,coordinates,order=1,mode='constant',cval=0.)
