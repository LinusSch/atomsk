#!/bin/bash

# Script to test as well as show example usage of the out_tric.f90 module
# for Atomsk.
#
# Created by Linus Schönström on 2020-08-12
# Last modified by Linus Schönström on 2020-08-13


# Changes to output directory and removes files from previous runs
cd output
rm *


# Translates the file NiSi2.cif to TRIC as well as two comparison points
atomsk ../input/NiSi2.cif \
	test_0 trc d12 sxyz 

	# This is testcase 0, the one that was used when developing the module
	# in 2020-08. Reference formats chosen:
	#
	# d12: fractional coordinates, 3 lengths + 3 angles for the cell
	#
	# sxyz: probable ground truth*, matrix representation of cell, can be
	# visualized with Vesta.
	#
	# *given its use when passing -v 4 (verbosity debug)



# The same, but sorting the atom list
atomsk ../input/NiSi2.cif \
	-sort z up -sort y up -sort x up \
	test_sorted trc d12 sxyz
	# This doesn't quite work as expected, despite sorting happening in the
	# defined order according to terminal output the written output files
	# have y sorted before z.



# Yanked from the examples already included with Atomsk
	## This script builds a carbon multi-wall nanotube (MWNT)
	## using atomsk. When creating a nanotube, its axis
	## is along Z and the nanotube is centered at (0,0).

	## Create (8,0) CNT
	atomsk --create nanotube 2.6 8 0 C cnt1.atsk

	## Create (16,0) CNT
	atomsk --create nanotube 2.6 16 0 C cnt2.atsk

	## Merge the two nanotubes, and repeat it 4 times along its axis (i.e. Z axis).
	## Note that the biggest NT is read first so that its cell encloses the two nanotubes
	atomsk --merge 2 cnt2.atsk cnt1.atsk mwnt.xsf sxyz trc -duplicate 1 1 4
		# This one should produce a warning and no output for TRIC, as the cell is not orthogonal

	# Then with -orthogonal-cell, this one should work with TRIC... well, not sure about the 
	# simulation code, but the output module should do it
	#
	# Note that this removes the negative coordinates in the above version so I need another test
	atomsk --merge 2 cnt2.atsk cnt1.atsk mwnt_orth.xsf sxyz trc -duplicate 1 1 4 -orthogonal-cell

	## Remove temporary files
	rm -f *.atsk
