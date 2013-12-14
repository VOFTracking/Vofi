# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([Vofi], [1.0], [simone.bna@unibo.it])

LIBS=-lm
USER_CFLAGS="$CFLAGS"

AC_PROG_CC
AC_PROG_CXX
AC_PROG_FC
AC_PROG_RANLIB
AC_PROG_INSTALL

AS_IF([test "$GCC" = yes],
  [AS_CASE([`$CC --version 2>&1 < /dev/null`],
    [*LLVM*], [opt=-O3],
    [*gcc*4.2* | *gcc*4.4.3*], [opt=-O0],
    [opt=-O3])
   CFLAGS=${USER_CFLAGS:-$opt -fomit-frame-pointer -ffast-math -Wall}],
  [CFLAGS=${USER_CFLAGS:--O}])
  
AC_C_CONST

AC_PATH_PROG(RM, rm, /bin/rm)

AC_PATH_PROG(AR, ar, no)
if [[ $AR = "no" ]] ; then
    AC_MSG_ERROR("Could not find ar - needed to create a library")
fi

AC_CONFIG_FILES([makefile])

AC_OUTPUT