from . import elements as pyblep_elements


def _from_madx_sequence(line, sequence, classes=pyblep_elements,
                        ignored_madtypes=[], exact_drift=False):

    if exact_drift:
        myDrift = classes.ExactDrift
    else:
        myDrift = classes.Drift
    seq = sequence

    elements = seq.elements
    ele_pos = seq.element_positions()

    old_pp = 0.
    i_drift = 0
    for ee, pp in zip(elements, ele_pos):
        if pp > old_pp:
            line.elements.append(myDrift(length=(pp-old_pp)))
            line.element_names.append('drift_%d' % i_drift)
            i_drift += 1

        eename = ee.name
        mad_etype = ee.base_type.name

        if ee.length > 0:
            raise ValueError(f"Sequence {seq} contains {eename} with length>0")


        if mad_etype in ['marker', 'monitor', 'hmonitor', 'vmonitor',
                         'rcollimator', 'placeholder', 'instrument', 'solenoid', 'drift']:
            newele = myDrift(length=ee.l)

        elif mad_etype == 'multipole':
            knl = ee.knl if hasattr(ee, 'knl') else [0]
            ksl = ee.ksl if hasattr(ee, 'ksl') else [0]
            newele = classes.Multipole(
                knl=knl,
                ksl=ksl,
                hxl=knl[0],
                hyl=0,
                length=ee.lrad)

        elif mad_etype == 'tkicker':
            hkick = [-ee.hkick] if hasattr(ee, 'hkick') else []
            vkick = [ee.vkick] if hasattr(ee, 'vkick') else []
            newele = classes.Multipole(knl=hkick, ksl=vkick,
                                       length=ee.lrad, hxl=0, hyl=0)

        elif mad_etype == 'vkicker':
            newele = classes.Multipole(knl=[], ksl=[ee.kick],
                                       length=ee.lrad, hxl=0, hyl=0)

        elif mad_etype == 'hkicker':
            newele = classes.Multipole(knl=[-ee.kick], ksl=[],
                                       length=ee.lrad, hxl=0, hyl=0)

        elif mad_etype == 'rfcavity':
            newele = classes.Cavity(voltage=ee.volt * 1e6,
                                    frequency=ee.freq * 1e6, lag=ee.lag * 360)

        elif mad_etype == 'beambeam':
            if ee.slot_id==6 or ee.slot_id==60:
                # BB interaction is 6D
                newele = classes.BeamBeam6D(
                    phi=0.,
                    alpha=0.,
                    x_bb_co=0.,
                    y_bb_co=0.,
                    charge_slices=[0.],
                    zeta_slices=[0.],
                    sigma_11=1.,
                    sigma_12=0.,
                    sigma_13=0.,
                    sigma_14=0.,
                    sigma_22=1.,
                    sigma_23=0.,
                    sigma_24=0.,
                    sigma_33=0.,
                    sigma_34=0.,
                    sigma_44=0.,
                    x_co=0.,
                    px_co=0.,
                    y_co=0.,
                    py_co=0.,
                    zeta_co=0.,
                    delta_co=0.,
                    d_x=0.,
                    d_px=0.,
                    d_y=0.,
                    d_py=0.,
                    d_zeta=0.,
                    d_delta=0.,
                )
            else:
                # BB interaction is 4D
                newele = classes.BeamBeam4D(
                    charge=0.,
                    sigma_x=1.,
                    sigma_y=1.,
                    beta_r=1.,
                    x_bb=0.,
                    y_bb=0.,
                    d_px=0.,
                    d_py=0.,
                )
        elif mad_etype in ignored_madtypes:
            pass

        else:
            raise ValueError('Not recognized')

        line.elements.append(newele)
        line.element_names.append(eename)

    other_info = {}

    return line, other_info
