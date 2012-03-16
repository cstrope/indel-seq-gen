#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
# in order to change the version of iSG, change the 2.0 below to the next version number.
AC_INIT(indel-seq-gen, 2.2.13, corystrope@gmail.com)
AC_CONFIG_SRCDIR([src/model.cpp])
AM_INIT_AUTOMAKE([-Wall -Werror])
AC_CONFIG_HEADER([config.h])

# Checks for programs.
AC_PROG_CC
gl_EARLY
AC_PROG_CXX
AC_PROG_CXXCPP

AC_LANG_PUSH(C++)

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lm':
AC_CHECK_LIB([m], [main])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([limits.h stdlib.h string.h unistd.h getopt.h ])
# For gnulib
gl_INIT

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_TYPE_SIZE_T
AC_C_INLINE
AC_C_BIGENDIAN

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([pow sqrt])

AX_CREATE_STDINT_H

AC_LANG_POP(C++)

AC_CONFIG_FILES([Makefile
		 lib/Makefile
		 data/Makefile
		doc/Makefile
         src/Makefile])

AC_OUTPUT
