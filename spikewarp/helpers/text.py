"""
ROUNDING HELPERS
Helper functions for rounding floats
"""

from math import log10, floor

def round_sig(x, sig=2):
	if (x == 0.0):
		return 0.0
	else:
		return round(x, sig-int(floor(log10(abs(x))))-1)
	
def round_to_1(x):
	return round_sig(x, sig=1)
def round_to_2(x):
	return round_sig(x, sig=2)
def round_to_3(x):
	return round_sig(x, sig=3)


def string_with_pval_lessthan_to_lower_limit(pvalue):

	if (pvalue < 0.001):
		return 'p < .001'
	else:
		"p = " + str(round_to_3(pvalue))



"""
PERCENTAGE & FRACTION HELPER
Helper function for creating string with percentage and corresponding fraction
"""

def percent_and_frac_string(numerator, denominator):
		
	percent = 0.0
	if (denominator > 0):
		percent = 100.0 * float(numerator) / float(denominator)
	percent_str = "{:.{}f}".format(percent, 1)
	full_str = percent_str + "\\% (" + str(numerator) + "/" + str(denominator) + ")" 

	return full_str


"""
LATEX HELPER
Helper function for writing latex tag to file
"""

def append_new_tag(text, tag_name, file_name):

	with open(file_name, "a") as tex_file:
		print("%<*" + tag_name + '>\n' + text + "\n%</" + tag_name + '>\n', file=tex_file)



