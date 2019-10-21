from scipy.special import wofz


class MathlibDefault(object):

    from numpy import sqrt, exp, sin, cos, abs, pi, tan

    @classmethod
    def wfun(cls, z_re, z_im):
        w = wofz(z_re + 1j * z_im)
        return w.real, w.imag
