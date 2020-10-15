all : lib Rfuncs

Rfuncs: src/dmatrix.h src/glasso.h 
	R CMD SHLIB src/glasso.c src/dmatrix.c -o lib/glasso.so
	R CMD SHLIB src/SFGen.c -o lib/SFGen.so

clean: 
	rm -f src/*.o src/*.so lib/*.so

lib: 
	mkdir lib
