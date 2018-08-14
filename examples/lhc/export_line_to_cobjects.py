import sixtracktools
import pysixtrack
from cobjects import CBuffer, CObject, CField

import pdb

def factorial( x ):
    if  not isinstance( x, ( int ) ):
        return 0
    return x > 0 and x * factorial( x - 1 ) or 1

class Drift( CObject ):
    _typeid =  2
    length  = CField( 0, 'real', default=0.0, alignment=8 )

class DriftExact( CObject ):
    _typeid = 3
    length  = CField( 0, 'real', default=0.0, alignment=8 )

class MultiPole( CObject ):
    _typeid = 4
    length  = None
    hxl     = None
    hyl     = None
    order   = None
    bal     = None

    def __init__( self, order = 1, cbuffer=None, **kwargs ):
        CObject.__init__( self, cbuffer=cbuffer, **kwargs )

        self.length = CField( 0, 'real',    default=0.0,  alignment=8 )
        self.hxl    = CField( 1, 'real',    default=0.0,  alignment=8 )
        self.hyl    = CField( 2, 'real',    default=0.0,  alignment=8 )
        self.order  = CField( 3, 'int64',   default=0,    alignment=8 )
        self.bal    = CField( 4, 'real',    default=0.0,
                                 length=( 2 * order + 2 ),
                                 pointer=True, alignment=8 )
        self.order  = order


class Cavity( CObject ):
    _typeid = 5
    voltage   = CField( 0, 'real', default=0.0,  alignment=8 )
    frequency = CField( 1, 'real', default=0.0,  alignment=8 )
    lag       = CField( 2, 'real', default=0.0,  alignment=8 )

class XYShift( CObject ):
    _typeid = 6
    dx      = CField( 0, 'real',   default=0.0,  alignment=8 )
    dy      = CField( 1, 'real',   default=0.0,  alignment=8 )

class SRotation( CObject ):
    _typeid = 7
    cos_z   = CField( 0, 'real',   default=1.0,  alignment=8 )
    sin_z   = CField( 1, 'real',   default=0.0,  alignment=8 )

class BeamBeam4D( CObject ):
    pass

class BeamBeam6D( CObject ):
    pass


if  __name__ == '__main__':
    from math import pi
    from math import sin, cos
    import numpy as np

    six = sixtracktools.SixTrackInput('.')
    line, rest, iconv = six.expand_struct(convert=pysixtrack.element_types)
    deg2rad = pi / 180.0

    b = CBuffer()

    for label, elem_type, elem in line:
        if  elem_type == 'Drift':
            e = Drift( cbuffer=b )
            e.length = elem.length

        elif elem_type == 'DriftExact':
            e = DriftExact( cbuffer=b )
            e.length = elem.length

        elif elem_type == 'Multipole':
            order = max( len( elem.knl ), len( elem.ksl ) )

            assert( order >= len( elem.knl ) )
            assert( order >= len( elem.ksl ) )

            knl = [ 0.0 for x in range( order ) ]
            knl[ :len( elem.knl ) ] = elem.knl

            ksl = [ 0.0 for x in range( order ) ]
            ksl[ :len( elem.ksl ) ] = elem.ksl

            bal = np.zeros( 2 * order )

            for i in range( order ):
                inv_factorial = float( 1.0 / factorial( i ) )
                bal[ 2 * i     ] = knl[ i ] * inv_factorial
                bal[ 2 * i + 1 ] = ksl[ i ] * inv_factorial

            if  order > 0: order -= 1

            e = MultiPole( order=order, cbuffer=b )
            e.length = elem.length
            e.hxl    = elem.hxl
            e.hyl    = elem.hyl
            e.order  = order
            e.bal    = bal


        elif elem_type == 'XYShift':
            e = XYShift( cbuffer=b )
            e.dx = elem.dx
            e.dy = elem.dy


        elif elem_type == 'SRotation':
            angle_rad = deg2rad * elem.angle
            cos_z     = cos( angle_rad )
            sin_z     = sin( angle_rad )
            e = SRotation( cbuffer=b )
            e.cos_z = cos_z
            e.sin_z = sin_z


        elif elem_type == 'Cavity':
            e = Cavity( cbuffer=b )
            e.voltage   = elem.voltage
            e.frequency = elem.frequency
            e.lag       = elem.lag

        else:
            print( elem_type )
            pdb.set_trace()

    b.to_file( './lhc_beam_elements.bin' )






