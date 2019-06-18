from .bases import Element

class Line(Element):
    _description = [
            ('elements','','List of elements',[]),
            ('element_names','', 'List of element names',[]))
    _extra = []

    @classmethod
    def from_sixinput(self,sixinput):
        pass

    @classmethod
    def from_madx_sequence(self,seq):
        pass

    def to_dict(self,keepextra=False):
        out={}
        out['elements']=[ el.to_dict(keepextra) for el in self.elements]
        out['element_names']=self.elements_names[:]
        return out

    def append_line(self, line):
        # Append the elements
        if type(line) is Line:
            # got a pysixtrack line
            self.elements += line.elements
        else:
            # got a different type of line (e.g. pybplep)
            for ee in line.elements:
                type_name = ee.__class__.__name__
                newele = element_types[type_name](**ee._asdict())
                self.elements.append(newele)

        # Append the names
        self.element_names += line.element_names

        assert(len(self.elements) == len(self.element_names))

        return self

    def track(self, p):
        for el in self.elements:
            el.track(p)

    def track_elem_by_elem(self, p):
        out = []
        for el in self.elements:
            out.append(p.copy())
            el.track(p)
        return out

    def find_closed_orbit(self, p0c, guess=[0.,0.,0.,0.,0.,0.],
            method='Nelder-Mead'):

        def _one_turn_map(coord):
            pcl = Particles(p0c=p0c)
            pcl.x = coord[0]
            pcl.px = coord[1]
            pcl.y = coord[2]
            pcl.py = coord[3]
            pcl.zeta = coord[4]
            pcl.delta = coord[5]

            self.track(pcl)
            coord_out = np.array(
                [pcl.x, pcl.px, pcl.y, pcl.py, pcl.sigma, pcl.delta])

            return coord_out

        def _CO_error(coord):
            return np.sum((_one_turn_map(coord) - coord)**2)

        if method == 'get_guess':
            res = type('', (), {})()
            res.x = guess
        else:
            import scipy.optimize as so
            res = so.minimize(_CO_error, np.array(
                guess), tol=1e-20, method=method)

        pcl = Particles(p0c=p0c)

        pcl.x = res.x[0]
        pcl.px = res.x[1]
        pcl.y = res.x[2]
        pcl.py = res.x[3]
        pcl.zeta = res.x[4]
        pcl.delta = res.x[5]

        return pcl

    def enable_beambeam(self):

        for ee in self.elements:
            if ee.__class__.__name__ in ['BeamBeam4D', 'BeamBeam6D']:
                ee.enabled = True

    def disable_beambeam(self):

        for ee in self.elements:
            if ee.__class__.__name__ in ['BeamBeam4D', 'BeamBeam6D']:
                ee.enabled = False


    def beambeam_store_closed_orbit_and_dipolar_kicks(self, particle_on_CO,
            separation_given_wrt_closed_orbit_4D=True,
            separation_given_wrt_closed_orbit_6D=True):

        self.disable_beambeam()
        closed_orbit = self.track_elem_by_elem(particle_on_CO)

        self.enable_beambeam()

        for ie, ee in enumerate(self.elements):
            ## to transfer to beambeam.py

            if ee.__class__.__name__ == 'BeamBeam4D':
                if separation_given_wrt_closed_orbit_4D:
                    ee.x_bb += closed_orbit[ie].x
                    ee.y_bb += closed_orbit[ie].y

                # Evaluate dipolar kick
                ptemp = closed_orbit[ie].copy()
                ptempin = ptemp.copy()

                ee.track(ptemp)

                ee.d_px = ptemp.px - ptempin.px
                ee.d_py = ptemp.py - ptempin.py

            elif ee.__class__.__name__ == 'BeamBeam6D':
                if not separation_given_wrt_closed_orbit_6D:
                    raise ValueError('Not implemented!')

                # Store closed orbit
                ee.x_co = closed_orbit[ie].x
                ee.px_co = closed_orbit[ie].px
                ee.y_co = closed_orbit[ie].y
                ee.py_co = closed_orbit[ie].py
                ee.zeta_co = closed_orbit[ie].zeta
                ee.delta_co = closed_orbit[ie].delta

                # Evaluate 6d kick on closed orbit
                ptemp = closed_orbit[ie].copy()
                ptempin = ptemp.copy()

                ee.track(ptemp)

                ee.d_x = ptemp.x - ptempin.x
                ee.d_px = ptemp.px - ptempin.px
                ee.d_y = ptemp.y - ptempin.y
                ee.d_py = ptemp.py - ptempin.py
                ee.d_zeta = ptemp.zeta - ptempin.zeta
                ee.d_delta = ptemp.delta - ptempin.delta







