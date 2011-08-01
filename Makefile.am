ACLOCAL_AMFLAGS=-I m4
if DEBUG
AM_CFLAGS = -g3 -O0 -DDEBUG
else
AM_CFLAGS = -O2
endif

if HAVE_DOXYGEN
    DOXYDIR = docs
endif 

SUBDIRS = src gen test $(DOXYDIR)
