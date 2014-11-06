dnl TVPI Library - A polyhedral domain for abstract interpretation.
dnl
dnl Copyright (C) 2003 Axel Simon
dnl
dnl This program is free software; you can redistribute it and/or modify
dnl it under the terms of the GNU General Public License as published by
dnl the Free Software Foundation; either version 2 of the License, or
dnl (at your option) any later version.
dnl
dnl This program is distributed in the hope that it will be useful,
dnl but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
dnl GNU General Public License for more details.
dnl
dnl You should have received a copy of the GNU General Public License
dnl along with this program; if not, write to the Free Software
dnl Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
dnl

dnl ######################################################################
dnl Process this file with autoconf to produce a configure script.
dnl ######################################################################

dnl Initialize and set the name of the (tar-)package. 
AC_INIT([Two Variables Abstract Domain],
	0.1.2, [A.Simon@kent.ac.uk], tvpi)

dnl Ensure that configure knows where the sources are.
AC_CONFIG_SRCDIR(Makefile.am)

dnl The Makefile is created by automake.
AM_INIT_AUTOMAKE(tvpi, 0.2.0)

dnl Create a header file with all the necessary flags.
AM_CONFIG_HEADER([config.h])

dnl Put a CVS revision stamp into the generated configure file.
AC_REVISION($Revision: 1.8 $)

dnl Check system type.
AC_CANONICAL_HOST

AC_PROG_RANLIB

dnl Use the C++ compiler to compile the tests.
AC_LANG(C++)

dnl Checks for programs.
AC_PROG_CC
AC_PROG_CXX
AC_PROG_CPP
AM_PROG_CC_C_O
AC_PROG_INSTALL
AC_PROG_LN_S
AC_CHECK_TOOL(AR,ar)

# Checks for libraries.
AC_LANG_PUSH(C)
AC_CHECK_LIB(gmp, __gmpz_init, ,[
  AC_MSG_ERROR([Cannot find the GNU multi-precision library.])])
AC_LANG_POP(C)

# Testing for main here is cheating a bit. At least the GNU linker
# falls over as soon as it can't find the library, hence the test
# fails when the library is not there.
AC_CHECK_LIB(gmpxx, main, ,[
  AC_MSG_ERROR([Cannot find the GNU multi-precision library.])])



# Checks for header files.
AC_LANG_PUSH(C)
AC_HEADER_STDC
AC_CHECK_HEADERS([stddef.h stdlib.h])
AC_LANG_POP(C)

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T
 
# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([sqrt])
                                                                                
# Substitue the results.
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
