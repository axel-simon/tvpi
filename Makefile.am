## TVPI Library - A polyhedral domain for abstract interpretation.
##
## Copyright (C) 2003 Axel Simon
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
##

############################################################
## Process this file with automake to produce Makefile.in ##
############################################################

## automake options:
## - build the object files where the sources are
AUTOMAKE_OPTIONS = 1.8 subdir-objects

# these subdirs contain sources
VPATH = $(srcdir):$(srcdir)/src:$(srcdir)/interface/C:$(srcdir)/tests

lib_LIBRARIES = libtvpi.a

libtvpi_a_SOURCES = src/memory.cpp \
		    src/interval.cpp \
		    src/planar.cpp \
		    src/polyhedron.cpp \
		    src/tvpi.cpp \
		    src/affine.cpp \
		    src/tvpi_c.cpp \
		    src/landmark.cpp


AM_CXXFLAGS = -Wall -Wno-long-long

# -fno-implicit-templates

## Headers that should be installed.
include_HEADERS = src/common.hh src/lincomponent.hh \
		  src/memory.hh src/interval.hh \
		  src/planar.hh src/polyhedron.hh \
	          src/tvpi.hh interface/C/tvpi_c.h \
		  src/affine.hh src/landmark.hh \
		  src/tvpiexception.hh

INCLUDES = -I$(srcdir)/interface/C -I$(srcdir)/src

## Tests.
LDADD = libtvpi.a

check_PROGRAMS = \
	tests/testInterval1 \
	tests/testInterval2 \
	tests/testInterval3 \
	tests/testInterval4 \
	tests/testInterval5 \
	tests/testInterval6 \
	tests/testInterval7 \
	tests/testInterval8 \
	tests/testConvergents1 \
	tests/testConvergents2 \
	tests/testDenseTvpiAugment1 \
	tests/testDenseTvpiAugment2 \
	tests/testDenseTvpiApproximateInequalities1 \
	tests/testDenseTvpiApproximateInequalities2 \
	tests/testDenseTvpiApproximateInequalities3 \
	tests/testDenseTvpiResultants2 \
	tests/testDenseTvpiRenameVariables1 \
	tests/testInequalityIsEntailed1 \
	tests/testInequalityIsEntailed2 \
	tests/testInequalityIsEntailed3 \
	tests/testInequalityCalculateCut1 \
	tests/testInequalityCalculateCut2 \
	tests/testInequalityCalculateCut3 \
	tests/testInequalityCalculateCut4 \
	tests/testInequalityCalculateCut5 \
	tests/testInequalityCalculateCut6 \
	tests/testInequalityCalculateCut7 \
	tests/testInequalityCalculateCut8 \
	tests/testInequalityCalcDirection1 \
	tests/testInequalityLessThanPi1 \
	tests/testPolyhedronWiden1 \
	tests/testPolyhedronWiden2 \
	tests/testPolyhedronWiden3 \
	tests/testPolyhedronHull1 \
	tests/testPolyhedronHull2 \
	tests/testPolyhedronHull3 \
	tests/testPolyhedronHull4 \
	tests/testPolyhedronHull5 \
	tests/testPolyhedronHull6 \
	tests/testPolyhedronHull7 \
	tests/testPolyhedronHull8 \
	tests/testPolyhedronHull9 \
	tests/testPolyhedronHull10 \
	tests/testPolyhedronAddInequalitySet1 \
	tests/testPolyhedronAddInequalitySet2 \
	tests/testPolyhedronAddInequalitySet3 \
	tests/testPolyhedronAddInequalitySet4 \
	tests/testPolyhedronAddInequalitySet5 \
	tests/testPolyhedronAddInequalitySet6 \
	tests/testPolyhedronAddInequalitySet7 \
	tests/testPolyhedronAddInequalitySet9 \
	tests/testPolyhedronPropagateBounds1 \
	tests/testDomainRemoveVariable2 \
	tests/testDomainUpdate1 \
	tests/testDomainEntails1 \
	tests/testDomainInequality1 \
	tests/testDomainInequality2 \
	tests/testDomainInequality3 \
	tests/testDomainJoin1 \
	tests/testDomainJoin2 \
	tests/testLandmarkTable1

#	tests/testDomainRemoveVariable1 
#	tests/testDenseTvpiApproximateInequalities3
#	tests/testDenseTvpiResultants1
#	tests/testDenseTvpiResultants3 
#	tests/testInequalityIsEntailed4 
#	tests/testPolyhedronSwap1 
#	tests/testPolyhedronIsSatisfiable1 
#	tests/testPolyhedronIsSatisfiable2
#	tests/testPolyhedronEntails1 
#	tests/testPolyhedronEntails2 
#	tests/testPolyhedronEntails3 
#	tests/testPolyhedronAddInequalitySet8 
#	tests/testPolyhedronEnforceBounds1 

TESTS = $(check_PROGRAMS)

tests_testInterval1_SOURCES = \
	tests/testInterval1.cpp
tests_testInterval2_SOURCES = \
	tests/testInterval2.cpp
tests_testInterval3_SOURCES = \
	tests/testInterval3.cpp
tests_testInterval4_SOURCES = \
	tests/testInterval4.cpp
tests_testInterval5_SOURCES = \
	tests/testInterval5.cpp
tests_testInterval6_SOURCES = \
	tests/testInterval6.cpp
tests_testInterval7_SOURCES = \
	tests/testInterval7.cpp
tests_testInterval8_SOURCES = \
	tests/testInterval8.cpp
tests_testConvergents1_SOURCES = \
	tests/testConvergents1.cpp
tests_testConvergents2_SOURCES = \
	tests/testConvergents2.cpp
tests_testDenseTvpiAugment1_SOURCES = \
	tests/testDenseTvpiAugment1.cpp
tests_testDenseTvpiAugment2_SOURCES = \
	tests/testDenseTvpiAugment2.cpp
tests_testDenseTvpiApproximateInequalities1_SOURCES = \
	tests/testDenseTvpiApproximateInequalities1.cpp
tests_testDenseTvpiApproximateInequalities2_SOURCES = \
	tests/testDenseTvpiApproximateInequalities2.cpp
tests_testDenseTvpiApproximateInequalities3_SOURCES = \
	tests/testDenseTvpiApproximateInequalities3.cpp
#tests_testDenseTvpiResultants1_SOURCES = \
#	tests/testDenseTvpiResultants1.cpp
tests_testDenseTvpiResultants2_SOURCES = \
	tests/testDenseTvpiResultants2.cpp
#tests_testDenseTvpiResultants3_SOURCES = \
#	tests/testDenseTvpiResultants3.cpp
tests_testDenseTvpiRenameVariables1_SOURCES = \
	tests/testDenseTvpiRenameVariables1.cpp
tests_testInequalityIsEntailed1_SOURCES = \
	tests/testInequalityIsEntailed1.cpp
tests_testInequalityIsEntailed2_SOURCES = \
	tests/testInequalityIsEntailed2.cpp
tests_testInequalityIsEntailed3_SOURCES = \
	tests/testInequalityIsEntailed3.cpp
#tests_testInequalityIsEntailed4_SOURCES = \
#	tests/testInequalityIsEntailed4.cpp
tests_testInequalityCalculateCut1_SOURCES = \
	tests/testInequalityCalculateCut1.cpp
tests_testInequalityCalculateCut2_SOURCES = \
	tests/testInequalityCalculateCut2.cpp
tests_testInequalityCalculateCut3_SOURCES = \
	tests/testInequalityCalculateCut3.cpp
tests_testInequalityCalculateCut4_SOURCES = \
	tests/testInequalityCalculateCut4.cpp
tests_testInequalityCalculateCut5_SOURCES = \
	tests/testInequalityCalculateCut5.cpp
tests_testInequalityCalculateCut6_SOURCES = \
	tests/testInequalityCalculateCut6.cpp
tests_testInequalityCalculateCut7_SOURCES = \
	tests/testInequalityCalculateCut7.cpp
tests_testInequalityCalculateCut8_SOURCES = \
	tests/testInequalityCalculateCut8.cpp
tests_testInequalityCalcDirection1_SOURCES = \
	tests/testInequalityCalcDirection1.cpp
tests_testInequalityLessThanPi1_SOURCES = \
	tests/testInequalityLessThanPi1.cpp
#tests_testPolyhedronSwap1_SOURCES = \
#	tests/testPolyhedronSwap1.cpp
#tests_testPolyhedronIsSatisfiable1_SOURCES = \
#	tests/testPolyhedronIsSatisfiable1.cpp
#tests_testPolyhedronIsSatisfiable2_SOURCES = \
#	tests/testPolyhedronIsSatisfiable2.cpp
#tests_testPolyhedronEntails1_SOURCES = \
#	tests/testPolyhedronEntails1.cpp
#tests_testPolyhedronEntails2_SOURCES = \
#	tests/testPolyhedronEntails2.cpp
#tests_testPolyhedronEntails3_SOURCES = \
#	tests/testPolyhedronEntails3.cpp
tests_testPolyhedronWiden1_SOURCES = \
	tests/testPolyhedronWiden1.cpp
tests_testPolyhedronWiden2_SOURCES = \
	tests/testPolyhedronWiden2.cpp
tests_testPolyhedronWiden3_SOURCES = \
	tests/testPolyhedronWiden3.cpp
tests_testPolyhedronHull1_SOURCES = \
	tests/testPolyhedronHull1.cpp
tests_testPolyhedronHull2_SOURCES = \
	tests/testPolyhedronHull2.cpp
tests_testPolyhedronHull3_SOURCES = \
	tests/testPolyhedronHull3.cpp
tests_testPolyhedronHull4_SOURCES = \
	tests/testPolyhedronHull4.cpp
tests_testPolyhedronHull5_SOURCES = \
	tests/testPolyhedronHull5.cpp
tests_testPolyhedronHull6_SOURCES = \
	tests/testPolyhedronHull6.cpp
tests_testPolyhedronHull7_SOURCES = \
	tests/testPolyhedronHull7.cpp
tests_testPolyhedronHull8_SOURCES = \
	tests/testPolyhedronHull8.cpp
tests_testPolyhedronHull9_SOURCES = \
	tests/testPolyhedronHull9.cpp
tests_testPolyhedronHull10_SOURCES = \
	tests/testPolyhedronHull10.cpp
tests_testPolyhedronAddInequalitySet1_SOURCES = \
	tests/testPolyhedronAddInequalitySet1.cpp
tests_testPolyhedronAddInequalitySet2_SOURCES = \
	tests/testPolyhedronAddInequalitySet2.cpp
tests_testPolyhedronAddInequalitySet3_SOURCES = \
	tests/testPolyhedronAddInequalitySet3.cpp
tests_testPolyhedronAddInequalitySet4_SOURCES = \
	tests/testPolyhedronAddInequalitySet4.cpp
tests_testPolyhedronAddInequalitySet5_SOURCES = \
	tests/testPolyhedronAddInequalitySet5.cpp
tests_testPolyhedronAddInequalitySet6_SOURCES = \
	tests/testPolyhedronAddInequalitySet6.cpp
tests_testPolyhedronAddInequalitySet7_SOURCES = \
	tests/testPolyhedronAddInequalitySet7.cpp
#tests_testPolyhedronAddInequalitySet8_SOURCES = \
#	tests/testPolyhedronAddInequalitySet8.cpp
tests_testPolyhedronAddInequalitySet9_SOURCES = \
	tests/testPolyhedronAddInequalitySet9.cpp
tests_testPolyhedronPropagateBounds1_SOURCES = \
	tests/testPolyhedronPropagateBounds1.cpp
#tests_testPolyhedronEnforceBounds1_SOURCES = \
#	tests/testPolyhedronEnforceBounds1.cpp
#tests_testDomainRemoveVariable1_SOURCES = \
#	tests/testDomainRemoveVariable1.cpp
tests_testDomainRemoveVariable2_SOURCES = \
	tests/testDomainRemoveVariable2.cpp
tests_testDomainUpdate1_SOURCES = \
	tests/testDomainUpdate1.cpp
tests_testDomainEntails1_SOURCES = \
	tests/testDomainEntails1.cpp
tests_testDomainInequality1_SOURCES = \
	tests/testDomainInequality1.cpp
tests_testDomainInequality2_SOURCES = \
	tests/testDomainInequality2.cpp
tests_testDomainInequality3_SOURCES = \
	tests/testDomainInequality3.cpp
tests_testDomainJoin1_SOURCES = \
	tests/testDomainJoin1.cpp
tests_testDomainJoin2_SOURCES = \
	tests/testDomainJoin2.cpp
tests_testLandmarkTable1_SOURCES = \
	tests/testLandmarkTable1.cpp


