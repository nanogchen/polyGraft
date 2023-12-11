#!/bin/env python

import sys
sys.path.insert(0,"../polyGraft/")
from polymer import cgPolymer

if __name__ == '__main__':

	# generate a linear polymer chain
	linear = cgPolymer()
	linear.setChain(6, [1.0])

	# save data
	linear.toDATA(f"linear_N6.data")	
	
