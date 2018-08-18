import numpy as np

def get_elems_of_type(line, elemname):

	pyst_ele_type_list = [ele[1] for ele in line]
	indlist = np.where([ss==elemname for ss in pyst_ele_type_list])[0]
	ele_list = [line[ind][2] for ind in indlist]
	name_list = [line[ind][0] for ind in indlist]

	return indlist, name_list, ele_list