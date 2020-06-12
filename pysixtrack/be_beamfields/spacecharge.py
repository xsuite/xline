from pysixtrack.base_classes import Element
from .gaussian_fields import get_Ex_Ey_Gx_Gy_gauss


class SpaceChargeCoasting(Element):
    """Space charge for a coasting beam."""

    _description = [
        ("number_of_particles", "", "Number of particles in the beam", 0.0),
        ("circumference", "m", "Machine circumference", 0.0),
        ("sigma_x", "m", "Horizontal size of the beam (r.m.s.)", 1.0),
        ("sigma_y", "m", "Vertical size of the beam (r.m.s.)", 1.0),
        ("length", "m", "Integration length of space charge kick", 0.0),
        ("x_co", "m", "Horizontal closed orbit offset", 0.0),
        ("y_co", "m", "Vertical closed orbit offset", 0.0),
    ]
    _extra = [
        ("min_sigma_diff", "m", "Threshold to detect round beam", 1e-10),
        ("enabled", "", "Switch to disable space charge effect", True),
    ]

    def track(self, p):
        if self.enabled:
            length = self.length
            sigma_x = self.sigma_x
            sigma_y = self.sigma_y

            charge = p.q0 * p.echarge
            x = p.x - self.x_co
            px = p.px
            y = p.y - self.y_co
            py = p.py

            chi = p.chi

            beta = p.beta0 / p.rvv
            p0c = p.p0c * p.echarge

            Ex, Ey = get_Ex_Ey_Gx_Gy_gauss(
                x,
                y,
                sigma_x,
                sigma_y,
                min_sigma_diff=self.min_sigma_diff,
                skip_Gs=True,
                mathlib=p._m,
            )

            fact_kick = (
                chi
                * self.number_of_particles
                / self.circumference
                * (charge * p.qratio)
                * charge
                * (1 - p.beta0 * beta)
                / (p0c * beta)
                * length
            )

            px += fact_kick * Ex
            py += fact_kick * Ey

            p.px = px
            p.py = py


class SpaceChargeQGaussianProfile(Element):
    """Space charge for a bunched beam with generalised
    Gaussian profile.
    """

    _description = [
        ("number_of_particles", "", "Number of particles in the bunch", 0.0),
        ("bunchlength_rms", "m", "Length of the bunch (r.m.s.)", 1.0),
        ("sigma_x", "m", "Horizontal size of the beam (r.m.s.)", 1.0),
        ("sigma_y", "m", "Vertical size of the beam (r.m.s.)", 1.0),
        ("length", "m", "Integration length of space charge kick", 0.0),
        ("x_co", "m", "Horizontal closed orbit offset", 0.0),
        ("y_co", "m", "Vertical closed orbit offset", 0.0),
    ]
    _extra = [
        ("min_sigma_diff", "m", "Threshold to detect round beam", 1e-10),
        ("enabled", "", "Switch to disable space charge effect", True),
        ("q_parameter", "", "q parameter of generalised Gaussian distribution (q=1 for standard Gaussian)", 1.0),
        ("b_parameter", "", "b parameter of generalised Gaussian distribution (b=1 for standard Gaussian)", 1.0),
    ]

    def track(self, p):
        if self.enabled:
            pi = p._m.pi
            exp = p._m.exp
            sqrt = p._m.sqrt
            bunchlength_rms = self.bunchlength_rms
            length = self.length
            sigma_x = self.sigma_x
            sigma_y = self.sigma_y

            charge = p.q0 * p.echarge
            x = p.x - self.x_co
            px = p.px
            y = p.y - self.y_co
            py = p.py
            sigma = p.sigma

            chi = p.chi

            beta = p.beta0 / p.rvv
            p0c = p.p0c * p.echarge

            Ex, Ey = get_Ex_Ey_Gx_Gy_gauss(
                x,
                y,
                sigma_x,
                sigma_y,
                min_sigma_diff=self.min_sigma_diff,
                skip_Gs=True,
                mathlib=p._m,
            )

            fact_kick = (
                chi
                * (charge * p.qratio)
                * charge
                * (1 - p.beta0 * beta)
                / (p0c * beta)
                * length
            )

            fact_kick *= (
                self.number_of_particles
                / (bunchlength_rms * sqrt(2 * pi))
                * exp(
                    -0.5
                    * (sigma / bunchlength_rms)
                    * (sigma / bunchlength_rms)
                )
            )

            px += fact_kick * Ex
            py += fact_kick * Ey

            p.px = px
            p.py = py


class SpaceChargeInterpolatedProfile(Element):
    """Space charge for a bunched beam with discretised profile."""

    _description = [
        ("number_of_particles", "", "Number of particles in the bunch", 0.0),
        ("line_density_profile", "1/m", "Discretised list of density values with integral normalised to 1", [1.0]),
        ("dz", "m", "Unit distance between profile points", 0.0),
        ("z0", "m", "Start position of line density profile", 0.0),
        ("sigma_x", "m", "Horizontal size of the beam (r.m.s.)", 1.0),
        ("sigma_y", "m", "Vertical size of the beam (r.m.s.)", 1.0),
        ("length", "m", "Integration length of space charge kick", 0.0),
        ("x_co", "m", "Horizontal closed orbit offset", 0.0),
        ("y_co", "m", "Vertical closed orbit offset", 0.0),
    ]
    _extra = [
        ("min_sigma_diff", "m", "Threshold to detect round beam", 1e-10),
        ("enabled", "", "Switch to disable space charge effect", True),
    ]

    def track(self, p):
        if self.enabled:
            pi = p._m.pi
            exp = p._m.exp
            sqrt = p._m.sqrt
            bunchlength_rms = self.bunchlength_rms
            length = self.length
            sigma_x = self.sigma_x
            sigma_y = self.sigma_y

            charge = p.q0 * p.echarge
            x = p.x - self.x_co
            px = p.px
            y = p.y - self.y_co
            py = p.py
            sigma = p.sigma

            chi = p.chi

            beta = p.beta0 / p.rvv
            p0c = p.p0c * p.echarge

            Ex, Ey = get_Ex_Ey_Gx_Gy_gauss(
                x,
                y,
                sigma_x,
                sigma_y,
                min_sigma_diff=self.min_sigma_diff,
                skip_Gs=True,
                mathlib=p._m,
            )

            fact_kick = (
                chi
                * (charge * p.qratio)
                * charge
                * (1 - p.beta0 * beta)
                / (p0c * beta)
                * length
            )

            fact_kick *= (
                self.number_of_particles
                / (bunchlength_rms * sqrt(2 * pi))
                * exp(
                    -0.5
                    * (sigma / bunchlength_rms)
                    * (sigma / bunchlength_rms)
                )
            )

            px += fact_kick * Ex
            py += fact_kick * Ey

            p.px = px
            p.py = py
