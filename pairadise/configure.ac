AC_INIT(pairadise2, 1.0.0)

# Allow the user to specify the path for Python
AC_ARG_WITH(python, [ --with-python       specify the location of python ])

# Allow the user to specify where STAR is
AC_ARG_WITH(STAR, [ --with-STAR       specify the location of STAR ])

# Check for python -- after this, PATH_TO_PYTHON will either contain a valid
# path to a location of the python executable, or (in the case where python can
# be found in the $PATH) will be empty.
with_python=`echo "$with_python" | sed "s|~|${HOME}|g"` # expand tilde
if test "$with_python" == ""
then
  with_python="$PATH"
  AC_SUBST(PATH_TO_PYTHON, "")
else
  AC_SUBST(PATH_TO_PYTHON, $with_python)
fi
AC_CHECK_PROG(PYTHON_CHECK,python,yes,no,$with_python)
if test x"$PYTHON_CHECK" != x"yes" ; then
  AC_MSG_ERROR([Could not find python, configuration failed])
fi

# Check for STAR -- after this, PATH_TO_STAR will either contain a valid
# path to a location containing the STAR executable, or (in the case where it
# can be found in the $PATH) will be empty.
with_STAR=`echo "$with_STAR" | sed "s|~|${HOME}|g"` # expand tilde
if test "$with_STAR" == ""
then
  with_STAR="$PATH"
  AC_SUBST(PATH_TO_STAR, "")
else
  AC_SUBST(PATH_TO_STAR, $with_STAR)
fi
AC_CHECK_PROG(STAR_CHECK,STAR,yes,no,$with_STAR)
if test x"$STAR_CHECK" != x"yes" ; then
  AC_MSG_ERROR([Could not find STAR, configuration failed])
fi

# set language to cpp
AC_PROG_CXX
AC_LANG(C++)

# Check for python
AC_CHECK_PROG(PYTHON_CHECK,$with_python,yes)
if test x"$PYTHON_CHECK" != x"yes" ; then
  AC_MSG_ERROR([Could not find the python executable, configuration failed])
fi

# Check for STAR
AC_CHECK_PROG(STAR_CHECK,$with_STAR,yes)
if test x"$STAR_CHECK" != x"yes" ; then
  AC_MSG_ERROR([Could not find the STAR executable, configuration failed])
fi

# Generate files from templates
AC_OUTPUT(Makefile)
