check_PROGRAMS = test_gf2matrix \
    test_linear_ops \
    test_nonlinear_ops \
    test_aes\
    test_wb_aes \
    test_wb_aes_inv
test_libs = $(top_builddir)/src/libwb_aes.a \
    $(top_builddir)/gen/libwb_aes_gen.a \
    $(top_builddir)/src/libwb_aes.a
test_includes = -I$(top_srcdir)/gen -I$(top_srcdir)/src

test_gf2matrix_SOURCES = test_gf2matrix.c test_util.c
test_gf2matrix_CPPFLAGS = ${test_includes}
test_gf2matrix_LDADD = ${test_libs}

test_linear_ops_SOURCES = test_linear_ops.c test_util.c
test_linear_ops_CPPFLAGS = ${test_includes}
test_linear_ops_LDADD = ${test_libs}

test_nonlinear_ops_SOURCES = test_nonlinear_ops.c test_util.c
test_nonlinear_ops_CPPFLAGS = ${test_includes}
test_nonlinear_ops_LDADD = ${test_libs}

test_aes_SOURCES = test_aes.c test_util.c
test_aes_CPPFLAGS = ${test_includes}
test_aes_LDADD = ${test_libs}

test_wb_aes_SOURCES = test_wb_aes.c test_util.c
test_wb_aes_CPPFLAGS = ${test_includes}
test_wb_aes_LDADD = ${test_libs}

test_wb_aes_inv_SOURCES = test_wb_aes_inv.c test_util.c
test_wb_aes_inv_CPPFLAGS = ${test_includes}
test_wb_aes_inv_LDADD = ${test_libs}

TESTS = ${check_PROGRAMS}

EXTRA_DIST = test.h test_util.h