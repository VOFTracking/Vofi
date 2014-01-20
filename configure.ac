# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([Vofi], [1.0], [simone.bna@unibo.it])

AC_MSG_RESULT(---------------------------------------------)
AC_MSG_RESULT(----------- Configuring vofi ----------------)
AC_MSG_RESULT(---------------------------------------------)

AC_PROG_CC
AC_PROG_CXX
AC_PROG_FC
AC_PROG_RANLIB
AC_PROG_INSTALL

AC_C_CONST
AC_C_INLINE

AS_IF([test "x${ac_cv_env_CFLAGS_set}" = "x"],
      [CFLAGS=""])
AS_IF([test "x${ac_cv_env_FCFLAGS_set}" = "x"],
      [FCFLAGS=""])
AS_IF([test "x${ac_cv_env_CXXFLAGS_set}" = "x"],
      [CXXFLAGS=""])

AC_MSG_CHECKING([C compiler flags...])
AS_IF([test "x${ac_cv_env_CFLAGS_set}" = "x"],[
  AS_IF([test "$GCC" = "yes"],[
    CFLAGS="${CFLAGS} -O3 -fomit-frame-pointer -ffast-math -Wall"
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
     
AC_PATH_PROG(RM, rm, /bin/rm)

AC_PATH_PROG(AR, ar, no)
if [[ $AR = "no" ]] ; then
    AC_MSG_ERROR("Could not find ar - needed to create a library")
fi

LIBS=-lm

AC_CONFIG_FILES([makefile])

AC_OUTPUT

AC_MSG_RESULT(---------------------------------------------)
AC_MSG_RESULT(--------- Done Configuring vofi -------------)
AC_MSG_RESULT(---------------------------------------------)