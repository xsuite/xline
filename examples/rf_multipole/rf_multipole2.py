import xline


rfmult = {"freq": 100, "knl": [0.1, 0.2], "ksl": [0.3, 0.4]}
out = xline.mad_benchmark("rfmultipole", rfmult, x=0.001, y=0.0012, t=0.1)
mad, line, p_mad, p_six = out

mult = {"knl": [0.1, 0.2], "ksl": [0.3, 0.4]}
out = xline.mad_benchmark("multipole", mult, x=0, pt=0.0)
mad, line, p_mad, p_six = out
