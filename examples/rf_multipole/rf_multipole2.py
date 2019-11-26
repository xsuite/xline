import pysixtrack


rfmult = {"freq": 100, "knl": [0.1]}
out = pysixtrack.mad_benchmark("rfmultipole", rfmult, x=1, t=0.1)
mad, line, p_mad, p_six = out

mult = {"knl": [0.1]}
out = pysixtrack.mad_benchmark("multipole", mult, x=0, pt=0.1)
mad, line, p_mad, p_six = out
