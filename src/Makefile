CPPFLAGS += -g -Wall -fopenmp

GTEST_DIR = tst/googletest/googletest
GTEST_HEADERS = $(GTEST_DIR)/include/gtest/*.h \
				$(GTEST_DIR)/include/gtest/internal/*.h
GTEST_H = $(GTEST_DIR)/include/gtest/gtest.h

.PHONY: clean test

gtest-all.o : 
	@ echo "Compiling Google Test Framework (gtest-all)..."
	@ cd $(GTEST_DIR)/scripts/test && \
	make gtest-all.o && \
	mv gtest-all.o ../../../../../

gtest_main.o : 
	@ echo "Compiling Google Test Framework (gtest_main)"
	@ cd $(GTEST_DIR)/scripts/test && \
	make gtest_main.o && \
	mv gtest_main.o ../../../../../

kdtreetest.o : tst/kdtreetest.cpp kdtree.hpp $(GTEST_HEADERS)
	@ echo "Compiling KDTree Test..."
	@ g++ $(CPPFLAGS) -isystem $(GTEST_DIR)/include $^ -c

kdtreetest : kdtreetest.o kdtree.o gtest-all.o gtest_main.o
	g++ $(CPPFLAGS) -lpthread $^ -o $@

hydrologytest.o : tst/hydrologytest.cpp hydrology.hpp $(GTEST_HEADERS)
	@ echo "Compiling Hydrology Test..."
	@ g++ $(CPPFLAGS) -isystem $(GTEST_DIR)/include $^ -c

hydrologytest : hydrologytest.o hydrology.o kdtree.o gtest-all.o gtest_main.o
	g++ $(CPPFLAGS) -lpthread $^ -o $@

hydrofunctionstest.o : tst/hydrofunctionstest.cpp hydrologyFunctions.hpp 
	@ echo "Compiling Hydrology Fnctions Test..."
	@ g++ $(CPPFLAGS) -isystem $(GTEST_DIR)/include $^ -c

hydrofunctionstest : hydrofunctionstest.o hydrologyFunctions.o hydrology.o kdtree.o raster.o gtest-all.o gtest_main.o
	g++ $(CPPFLAGS) -lpthread -lopencv_imgproc $^ -o $@

test : kdtreetest hydrologytest hydrofunctionstest
	@ ./kdtreetest
	@ ./hydrologytest
	@ ./hydrofunctionstest

kdtree.o : kdtree.cpp kdtree.hpp point.hpp

hydrology.o : hydrology.cpp hydrology.hpp point.hpp

buildRivers.o : buildRivers.cpp hydrology.hpp kdtree.hpp raster.hpp

hydrologyParameters.o : hydrologyParameters.cpp raster.hpp hydrology.hpp

raster.o : raster.cpp raster.hpp

hydrologyFunctions.o : hydrologyFunctions.cpp hydrologyFunctions.hpp hydrology.hpp hydrologyParameters.hpp

buildRivers: hydrologyFunctions.o raster.o hydrologyParameters.o buildRivers.o hydrology.o kdtree.o
	g++ hydrologyFunctions.o buildRivers.o raster.o hydrologyParameters.o hydrology.o kdtree.o $(CPPFLAGS) -lopencv_core -lopencv_imgproc -lopencv_highgui -o $@

clean :
	@ rm -vf buildRivers kdtreetest hydrologytest hydrofunctionstest *.gch *.o
	@ cd $(GTEST_DIR)/scripts/test/ && make clean