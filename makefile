all:

clean:
	rm -Rf bin
	cmake -P scripts/clean.cmake
	${MAKE} -C src/lithium/doc clean
