all : lib opt

opt: src/dmatrix.h src/def.h src/opt.h
	R CMD SHLIB src/opt.cpp src/def.c src/dmatrix.c -o lib/opt.so
	R CMD SHLIB src/SFGen.c -o lib/SFGen.so

clean: 
	rm -f src/*.o src/*.so lib/*.so

lib: 
	mkdir lib
