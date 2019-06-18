# Read sixtrack input using sixtracktools
import sixtracktools
sixinput=sixtracktools.SixInput('sixtrack_input')

# Build a pysixtrack line from pyblep line
import pysixtrack
ps_line, other = pysixtrack.Line.from_sixinput(sixinput)

# # Build a pysixtracklib line from pyblep line
# import pysixtracklib
# pslib_line = pysixtracklib.Elements()
# pslib_line.append_line(line)
# pslib_line.BeamMonitor(num_stores=1)
# 
# # Build a pysixtrack particle
# ps_part = pysixtrack.Particles(p0c=7000e9)
# ps_part.x = 1e-3
# ps_part.px = 2e-4
# 
# # Build a pysixtracklib particle
# pslib_part_set = pysixtracklib.ParticlesSet()
# pslib_part = pslib_part_set.Particles(num_particles=1)
# ps_part.partid = 0
# ps_part.state = 1
# ps_part.elemid = 0
# ps_part.turn = 0
# pslib_part.from_pysixtrack(ps_part, particle_index=0)
# 
# # Track with pysisxtrack
# ps_line.track(ps_part)
# 
# # Track with sixtracklib
# job = pysixtracklib.TrackJob(pslib_line, pslib_part_set)
# job.track(until_turn=1)
# job.collect()
# 
# # Compare
# print('pysixtrack x: ', ps_part.x)
# print('sixtracklib x:',job.output.particles[0].x[0])
