#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_INIT([Vofi], [1.0], [simone.bna@unibo.it])
AC_PREREQ([2.69])

AC_MSG_RESULT(---------------------------------------------)
AC_MSG_RESULT(----------- Configuring vofi ----------------)
AC_MSG_RESULT(---------------------------------------------)

AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_AUX_DIR([.])

AC_CANONICAL_BUILD
AC_CANONICAL_HOST
AC_CANONICAL_TARGET

# Makefiles 
AC_CONFIG_FILES([makefile 
                 src/makefile 
                 demo_src/makefile
                 demo_src/C/2D/Ellipse/makefile
                 demo_src/C/2D/Gaussian/makefile
                 demo_src/C/2D/Rectangle/makefile
                 demo_src/C/2D/Sine_line/makefile
                 demo_src/C/3D/Cap1/makefile
                 demo_src/C/3D/Cap2/makefile
                 demo_src/C/3D/Cap3/makefile
                 demo_src/C/3D/Sine_surface/makefile
                 demo_src/C/3D/Sphere/makefile
                 demo_src/CPP/2D/Ellipse/makefile
                 demo_src/CPP/2D/Gaussian/makefile
                 demo_src/CPP/2D/Rectangle/makefile
                 demo_src/CPP/2D/Sine_line/makefile
                 demo_src/CPP/3D/Cap1/makefile
                 demo_src/CPP/3D/Cap2/makefile
                 demo_src/CPP/3D/Cap3/makefile
                 demo_src/CPP/3D/Sine_surface/makefile
                 demo_src/CPP/3D/Sphere/makefile
                 demo_src/Fortran/2D/Ellipse/makefile
                 demo_src/Fortran/2D/Gaussian/makefile
                 demo_src/Fortran/2D/Rectangle/makefile
                 demo_src/Fortran/2D/Sine_line/makefile
                 demo_src/Fortran/3D/Cap1/makefile
                 demo_src/Fortran/3D/Cap2/makefile
                 demo_src/Fortran/3D/Cap3/makefile
                 demo_src/Fortran/3D/Sine_surface/makefile
                 demo_src/Fortran/3D/Sphere/makefile
                ])

# require automake 
AM_INIT_AUTOMAKE

# use silent rules - automake 1.11
# AM_SILENT_RULES(yes)

# by default we want 'maintainer mode', which means automake will generate rules to recreate e.g. Makefile.in from
# Makefile.am. But this adds support for --disable-maintainer-mode, which may be useful for forcing configure to
# work on systems with no autotools
AM_MAINTAINER_MODE([disable])

# Checks for header files.
AC_HEADER_STDC
#AC_CHECK_HEADERS([stdlib.h])

# Checks for typedefs, structures, and compiler characteristics.

# Checks for libraries.
AC_CHECK_LIB([m],[cos])

# Checks for library functions.
AC_CHECK_FUNCS([sqrt])

# --------------------------------------------------------------
# Do compilation tests using the compiler, preprocessor, and file extensions
# of C.  This uses the .C extension for test programs.
# --------------------------------------------------------------
AC_LANG([C])

# Checks for programs.
AC_PROG_CC
AC_PROG_CXX
AC_PROG_FC
#AM_PROG_AR
AC_PATH_PROG(AR, ar, /bin/ar)
AC_PATH_PROG(RM, rm, /bin/rm)
RM="$RM -f"
AC_PROG_INSTALL
AC_C_CONST
AC_C_INLINE

# Compiler flags
# by default CFLAGS and CXXFLAGS are set to '-g -O2' on systems that support them.
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
#AC_SUBST([LIBTOOL_DEPS])
#AC_PROG_LIBTOOL
#AM_PROG_LIBTOOL

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

# Create output files. Also configures any subpackages
AC_OUTPUT

AC_MSG_RESULT()
AC_MSG_RESULT(-------------------------------------------------------------------------------)
AC_MSG_RESULT(Configure complete now type 'make' and then 'make install')
AC_MSG_RESULT()

AC_MSG_RESULT(---------------------------------------------)
AC_MSG_RESULT(--------- Done Configuring vofi -------------)
AC_MSG_RESULT(---------------------------------------------)
