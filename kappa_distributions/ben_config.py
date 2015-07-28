# Ben's frequently used constants and functions. To use in ipython, execute:
#
# import matplotlib
# try:
#   ben
# except NameError:
#   import ben_config as ben
# else:
#   ben = reload(ben)
# ben.latexify(matplotlib)

print 'ben_config loaded'
print 'Available constants:'
print '  CL, QE, ME, MP, MN, HPL, HBAR, KBOL, GNEWT, SIG, AR, THOMSON,'
print '  JY, PC, AU, YEAR, DAY, HOUR, MSOLAR, RSOLAR, LSOLAR'
print 'Available functions:'
print '  latexify(mpl)'
print '  planck(nu, Thetae) '
print '  TwoPointCorrelationFunction_2D(dataslice[Ny][Nx], Ny, Nx)'
import math
import numpy as np

CL      = 2.99792458e10
QE      = 4.80320680e-10
ME      = 9.1093826e-28
MP      = 1.67262171e-24
MN      = 1.67492728e-24
HPL     = 6.6260693e-27
HBAR    = 1.0545717e-27
KBOL    = 1.3806505e-16
GNEWT   = 6.6742e-8
SIG     = 5.670400e-5
AR      = 7.5657e-15
THOMSON = 0.665245873e-24
JY      = 1.e-23
PC      = 3.085678e18
AU      = 1.49597870691e13
YEAR    = 31536000.
DAY     = 86400.
HOUR    = 3600.
MSOLAR  = 1.989e33
RSOLAR  = 6.96e10
LSOLAR  = 3.827e33

def latexify(mpl):
    fig_width  = 15
    fig_height = (5.**0.5 - 1.)/2.*fig_width # Golden ratio

    params = {'backend'            : 'ps',
              'text.latex.preamble': ['\usepackage{gensymb}'],
              'axes.labelsize'     : 22,
              'axes.titlesize'     : 22,
              'text.fontsize'      : 22,
              'legend.fontsize'    : 22,
              'xtick.labelsize'    : 22,
              'ytick.labelsize'    : 22,
              'text.usetex'        : True,
              'figure.figsize'     : [fig_width,fig_height],
              'font.family'        : 'serif'
             }

    mpl.rcParams.update(params)

def planck(nu, Thetae):
    return 2.*HPL*nu**3./(CL*CL)*1./(math.expm1(HPL*nu/(ME*CL*CL*Thetae)))

def TwoPointCorrelationFunction_2D(dataslice, Ny, Nx):
    # FFT
    dataslice_fft = np.fft.fft2(dataslice, s=(Ny,Nx))
    # Square
    dataslice_squared = dataslice_fft[:,:].real**2. + dataslice_fft[:,:].imag**2.
    # Inverse FFT
    dataslice_inverted = np.fft.ifft2(dataslice_squared, s=(Ny,Nx))
    
    dataslice_inverted_shifted = np.zeros([Ny, Nx])
    for i in xrange(Nx):
        for j in xrange(Ny):
            ik = i + Nx/2
            jk = j + Ny/2
            if ik >= Nx:
                ik = ik - Nx
            if jk >= Ny:
                jk = jk - Ny
            dataslice_inverted_shifted[jk,ik] = dataslice_inverted[j,i].real/(Nx*Ny)
            
    return dataslice_inverted_shifted
