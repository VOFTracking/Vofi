#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_INIT([Vofi], [1.0], [simone.bna@unibo.it])
AC_PREREQ([2.69])

AC_MSG_RESULT(---------------------------------------------)
AC_MSG_RESULT(----------- Configuring vofi ----------------)
AC_MSG_RESULT(---------------------------------------------)

AC_CONFIG_MACRO_DIR([m4])

# --------------------------------------------------------------
# Compute the canonical target system type variable "target", and its
# three individual parts: "target_cpu", "target_vendor", and "target_os"
# --------------------------------------------------------------
AC_CANONICAL_TARGET

# --------------------------------------------------------------
# makefiles generation 
# --------------------------------------------------------------
AC_CONFIG_FILES([makefile 
                 src/makefile 
                 demo_src/makefile
                 demo_src/C/makefile
                 demo_src/CPP/makefile
                 demo_src/Fortran/makefile
                ])

# --------------------------------------------------------------
# require automake
# --------------------------------------------------------------
AM_INIT_AUTOMAKE

# --------------------------------------------------------------
# do not use silent rules
# --------------------------------------------------------------
AM_SILENT_RULES(no)

# --------------------------------------------------------------
# by default we do not want 'maintainer mode', which may be useful for forcing configure to
# work on systems with no autotools
# --------------------------------------------------------------
AM_MAINTAINER_MODE([disable])

# --------------------------------------------------------------
# Checks for header files.
# --------------------------------------------------------------
AC_HEADER_STDC

# --------------------------------------------------------------
# Checks for typedefs, structures, and compiler characteristics.
# --------------------------------------------------------------

# --------------------------------------------------------------
# Checks for libraries.
# --------------------------------------------------------------
AC_CHECK_LIB([m],[cos])

# --------------------------------------------------------------
# Checks for library functions.
# --------------------------------------------------------------
AC_CHECK_FUNCS([sqrt])

# --------------------------------------------------------------
# Do compilation tests using the compiler, preprocessor, and file extensions
# of C, CPP and Fotran.  This uses the .C extension for test programs.
# --------------------------------------------------------------

AC_PROG_CC
AC_PROG_CXX
AC_PROG_FC
AC_C_CONST
AC_C_INLINE
AC_LANG([C])
AC_LANG([C++])
AC_LANG([Fortran])


# --------------------------------------------------------------
# Set AR, RM and INSTALL variables
# --------------------------------------------------------------
AC_PATH_PROG(AR, ar, /bin/ar)
AC_PATH_PROG(RM, rm, /bin/rm)
RM="$RM -f"
AC_PROG_INSTALL

# --------------------------------------------------------------
# C, Cpp and Fortran compiler flags
# --------------------------------------------------------------
AS_IF([test "x${ac_cv_env_CFLAGS_set}" = "x"],
      [CFLAGS=""])
AS_IF([test "x${ac_cv_env_FCFLAGS_set}" = "x"],
      [FCFLAGS=""])
AS_IF([test "x${ac_cv_env_CXXFLAGS_set}" = "x"],
      [CXXFLAGS=""])

AC_MSG_CHECKING([C compiler flags...])
AS_IF([test "x${ac_cv_env_CFLAGS_set}" = "x"],[
  AS_IF([test "$GCC" = "yes"],[
    CFLAGS="${CFLAGS} -O2 -fomit-frame-pointer -ffast-math -Wall"
  ])
  AC_MSG_RESULT([use default: $CFLAGS])
],[
  AC_MSG_RESULT([use user-supplied: $CFLAGS])
])      

AC_MSG_CHECKING([Fortran compiler flags...])
AS_IF([test "x${ac_cv_env_FCFLAGS_set}" = "x"],[
  AS_IF([test "$GCC" = "yes"],[
    FCFLAGS="${FCFLAGS} -O2 -Wall"
  ])
  AC_MSG_RESULT([use default: $FCFLAGS])
],[
  AC_MSG_RESULT([use user-supplied: $FCFLAGS])
]) 
      
AC_MSG_CHECKING([C++ compiler flags...])
AS_IF([test "x${ac_cv_env_CXXFLAGS_set}" = "x"],[
  AS_IF([test "$GCC" = "yes"],[
    CXXFLAGS="${CXXFLAGS} -O2 -Wall"
  ])
  AC_MSG_RESULT([use default: $CXXFLAGS])
],[
  AC_MSG_RESULT([use user-supplied: $CXXFLAGS])
])

#-----------------------------------------------------
# Initialize libtool.  By default, we will build
# only shared libraries on platforms that support them
#-----------------------------------------------------
LT_INIT([disable-static])

# --------------------------------------------------------------
# Portable symbolic links
# --------------------------------------------------------------
AC_PROG_LN_S
# -------------------------------------------------------------

# --------------------------------------------------------------
# parent directories -- portably and without race conditions
# --------------------------------------------------------------
AC_PROG_MKDIR_P
# -------------------------------------------------------------

# --------------------------------------------------------------
# nice sed that doesn't truncate output
# --------------------------------------------------------------
AC_PROG_SED
AC_SUBST(SED)
# --------------------------------------------------------------

# --------------------------------------------------------------
# Create output files. 
# --------------------------------------------------------------
AC_OUTPUT

AC_MSG_RESULT()
AC_MSG_RESULT(-------------------------------------------------------------------------------)
AC_MSG_RESULT(Configure complete now type 'make' and then 'make install')
AC_MSG_RESULT()

AC_MSG_RESULT(---------------------------------------------)
AC_MSG_RESULT(--------- Done Configuring vofi -------------)
AC_MSG_RESULT(---------------------------------------------)
