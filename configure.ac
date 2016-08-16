AC_PREREQ(2.57)
AC_INIT([rtfbs], 0.3.4, [mjhubisz@cornell.edu], [], [http://compgen.cshl.edu/rtfbs])
#AC_CONFIG_SRCDIR([RTFBS/src/rph_init.c])

: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
   AC_MSG_ERROR([Could not determine R_HOME.])
fi

AC_MSG_CHECKING(for CFLAGS)
#AC_PROG_CC
CC=`"${R_HOME}/bin/R" CMD config CC`
#RCFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
RCFLAGS="-g "
#CFLAGS="${RCFLAGS} ${CFLAGS}"
CFLAGS="${CFLAGS} -g " # -fopenmp "
AC_MSG_RESULT(${CFLAGS})
AC_PROG_CC
AC_MSG_RESULT(${CFLAGS})

# Checks for header files.
AC_FUNC_ALLOCA
AC_HEADER_STDC
AC_CHECK_HEADERS([stddef.h stdlib.h string.h strings.h sys/time.h unistd.h])

AC_CHECK_HEADER(pcre.h, , [AC_MSG_ERROR([Couldn not find pcre.h... download/install at www.pcre.org])])


# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_HEADER_STDBOOL
AC_TYPE_SIZE_T
AC_HEADER_TIME
AC_STRUCT_TM

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_MEMCMP
AC_FUNC_REALLOC
AC_FUNC_STRTOD
AC_FUNC_VPRINTF
AC_CHECK_FUNCS([atexit bzero floor gethostname gettimeofday memset pow re_comp regcomp rint sqrt strcasecmp strdup strstr strtol])


AC_SEARCH_LIBS([pcre_compile], [pcre], [],
 [AC_MSG_RESULT(no)
  AC_MSG_ERROR([ Need to install pcre ])])

OLD_CFLAGS=${CFLAGS}

# removed -g flag (CRAN rules)
#CFLAGS="${CFLAGS} -fno-strict-aliasing "
CFLAGS="-fno-strict-aliasing" 
AC_MSG_CHECKING(for -fno-strict-aliasing)
AC_COMPILE_IFELSE([AC_LANG_PROGRAM()], AC_MSG_RESULT(yes), [AC_MSG_RESULT(no)
		                           CFLAGS=${OLD_CFLAGS}])

AC_C_INLINE

AC_SUBST(CFLAGS)

AC_OUTPUT(src/Makevars)

cat <<EOF

$PACKAGE-$VERSION configuration summary:

    C preprocessor .................. : ${CPP}
    C compiler ...................... : ${CC}
    C compiler flags ................ : ${CFLAGS}
    Extra libraries ................. : ${LIBS}

EOF

dnl end configure.ac

