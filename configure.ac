#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.64])
AC_INIT([WB_AES],[1.0],[wbaes@hotmail.com])
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([config.h.in])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])

# Checks for programs.
AC_PROG_CC
AM_PROG_CC_C_O
AC_PROG_INSTALL
AC_PROG_RANLIB
AC_PROG_MKDIR_P

AC_CHECK_PROGS([DOXYGEN], [doxygen])
if test -z "$DOXYGEN";
   then AC_MSG_WARN([Doxygen not found - continuing without Doxygen support])
fi
AM_CONDITIONAL([HAVE_DOXYGEN],
    [test -n "$DOXYGEN"])
AM_COND_IF([HAVE_DOXYGEN], 
    [AC_CONFIG_FILES([docs/doxyfile])])

# Checks for header files.
AC_CHECK_HEADERS([memory.h stdint.h stdlib.h string.h])

# Checks for libraries
# Check for M4RI static library
have_m4ri=no
AC_SEARCH_LIBS([mzd_init], [m4ri], [have_m4ri=yes])
if test "${have_m4ri}" = "yes"; then
	AC_CHECK_HEADER([m4ri/m4ri.h], [], [have_m4ri=no], [-])
fi
if test "${have_m4ri}" = "no"; then
	AC_MSG_ERROR([M4RI library and include path are needed.	Install M4RI and rerun configure])
fi
# End of check for M4RI

AC_ARG_ENABLE(debug,
AS_HELP_STRING([--enable-debug],
               [enable debugging, default: no]),
[case "${enableval}" in
             yes) debug=true ;;
             no)  debug=false ;;
             *)   AC_MSG_ERROR([bad value ${enableval} for --enable-debug]) ;;
esac],
[debug=false])

AM_CONDITIONAL(DEBUG, test x"$debug" = x"true")

# Checks for typedefs, structures, and compiler characteristics
AC_TYPE_UINT32_T
AC_TYPE_UINT8_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([memset memcpy])

AC_CONFIG_FILES([Makefile
    src/Makefile
    gen/Makefile
    test/Makefile
    docs/Makefile])

AC_OUTPUT

