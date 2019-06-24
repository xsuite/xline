from cpymad.madx import Madx
import pysixtrack

mad = Madx()
mad.options.echo = False
mad.options.warn = False
mad.options.info = False

mad.call('mad/lhcwbb.seq')

line, other = pysixtrack.Line.from_madx_sequence(mad.sequence.lhcb1)

print(line)
