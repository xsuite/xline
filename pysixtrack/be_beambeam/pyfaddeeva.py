from scipy.special import wofz

def wfun(z_re, z_im):
    w = wofz(z_re+1j*z_im)
    return w.real, w.imag
