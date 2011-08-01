
all-local: doxyfile.stamp

doxyfile.stamp: doxyfile
	$(DOXYGEN) $<
	echo Timestamp > $@

CLEANFILES = doxyfile.stamp

clean-local:
	rm -rf html
