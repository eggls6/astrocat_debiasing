import numpy as np
import pandas as pd
import healpy as hp


__all__ = ['debiasRADec','icrf2radec','radec2icrf']


def debiasRADec(ra,dec,epoch,catalog,biasdf,J2000=2451545.0,nside=256):
    
    """Astrometric catalog bias correction following Eggl et al. (2020). 
    
    Parameters:
    -----------
    ra      ... [rad] Right Ascension
    dec     ... [rad] Declination
    epoch   ... [JD]  epoch of observation (Julian Date)
    catalog ... [string] MPC catalog identifier (one letter)
    biasdf  ... [pandas] pandas DataFrame containing debiasing information 
    J2000   ... [JD] J2k epoch
    nside   ... [int] nside resolution for healpix tesselation (nside = 2^N)  
    
    Dependencies:
    -------------
    numpy, healpy, pandas
    
    Returns:
    --------
    ra_deb  ... [rad] debiased Right Ascension at epoch
    dec_deb ... [rad] debiased declination at epoch
    """
    # arcseconds to degrees
    as2deg = 1/3600
    # milliarcseconds to degrees
    mas2deg = 1/3600/1000
    # numpy convenience functions
    pi = np.pi
    pix2 = pi*2
    piby2 = pi/2
    cos = np.cos
    # healpix function mapping RADEC to pixel number
    ang2pix = hp.ang2pix
    
    # find pixel from RADEC
    idx = ang2pix(nside, ra, dec, nest=False, lonlat=True)
    # find catalog data in pandas Data Frame
    colnames = [catalog+'_ra',catalog+'_dec',catalog+'_pm_ra',catalog+'_pm_dec']
    bias=biasdf[colnames].iloc[idx]
    
    # time from epoch in Julian years
    dt=(epoch-J2000)/365
    
    # bias correction
    ddec = (bias[colnames[1]]*as2deg+dt*bias[colnames[3]]*mas2deg)
    dec_deb = dec-ddec
    
    dra = (bias[colnames[0]]*as2deg+dt*bias[colnames[2]]*mas2deg)
    ra_deb = ra-dra/cos(dec)
    
    # Quadrant correction
    xyz=radec2icrf(ra_deb, dec_deb, deg=False)  
    ra_deb, dec_deb = icrf2radec(xyz[0], xyz[1], xyz[2], deg=False)
    
    return ra_deb, dec_deb 


def icrf2radec(x, y, z, deg=True):
    """Convert ICRF xyz to Right Ascension and Declination.
    Geometric states on unit sphere, no light travel time/aberration correction.
    Parameters:
    -----------
    x,y,z ... 3D vector of unit length (ICRF)
    deg ... True: angles in degrees, False: angles in radians
    Returns:
    --------
    ra ... Right Ascension [deg]
    dec ... Declination [deg]
    """

    norm=np.linalg.norm
    array=np.array
    arctan2=np.arctan2
    arcsin=np.arcsin
    rad2deg=np.rad2deg
    modulo=np.mod
    pix2=2.*np.pi

    pos=array([x,y,z])
    if(pos.ndim>1):
        r=norm(pos,axis=0)
    else:
        r=norm(pos)

    xu=x/r
    yu=y/r
    zu=z/r

    phi=arctan2(yu,xu)
    delta=arcsin(zu)

    if(deg):
        ra = modulo(rad2deg(phi)+360,360)
        dec = rad2deg(delta)
    else:
        ra = modulo(phi+pix2,pix2)
        dec = delta

    return ra, dec


def radec2icrf(ra, dec, deg=True):
    """Convert Right Ascension and Declination to ICRF xyz unit vector.
    Geometric states on unit sphere, no light travel time/aberration correction.
    Parameters:
    -----------
    ra ... Right Ascension [deg]
    dec ... Declination [deg]
    deg ... True: angles in degrees, False: angles in radians
    Returns:
    --------
    x,y,z ... 3D vector of unit length (ICRF)
    """

    deg2rad=np.deg2rad
    array=np.array
    cos=np.cos
    sin=np.sin

    if(deg):
        a = deg2rad(ra)
        d = deg2rad(dec)
    else:
        a = array(ra)
        d = array(dec)

    cosd = cos(d)
    x = cosd*cos(a)
    y = cosd*sin(a)
    z = sin(d)

    return array([x, y, z])
