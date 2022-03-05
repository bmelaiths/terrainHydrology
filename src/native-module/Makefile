CPPFLAGS += -g -Wall -fopenmp -isystem /usr/include/opencv4

BIN = ./bin
GEN = $(BIN)/gen

GTEST_DIR = tst/googletest/googletest
GTEST_HEADERS = $(GTEST_DIR)/include/gtest/*.h \
				$(GTEST_DIR)/include/gtest/internal/*.h
GTEST_H = $(GTEST_DIR)/include/gtest/gtest.h

.PHONY: clean test all

all: $(BIN)/buildRivers $(BIN)/terrainPrimitives

# gtest-all.o : 
# 	@ echo "Compiling Google Test Framework (gtest-all)..."
# 	@ cd $(GTEST_DIR)/scripts/test && \
# 	make gtest-all.o && \
# 	mv gtest-all.o ../../../../../

# gtest_main.o : 
# 	@ echo "Compiling Google Test Framework (gtest_main)"
# 	@ cd $(GTEST_DIR)/scripts/test && \
# 	make gtest_main.o && \
# 	mv gtest_main.o ../../../../../

# testsuite.o : tst/testsuite.cpp kdtree.hpp raster.hpp forest.hpp $(GTEST_HEADERS)
# 	@ echo "Compiling Tests..."
# 	@ g++ $(CPPFLAGS) -isystem $(GTEST_DIR)/include $^ -c

# testsuite : testsuite.o gtest-all.o gtest_main.o terrainHoneycomb.o terrainElevation.o hydrologyFunctions.o hydrologyParameters.o hydrology.o floatEndian.o shore.o
# 	g++ $(CPPFLAGS) $^ -lpthread -lopencv_core -lopencv_imgproc -lgeos_c -o $@

# test : testsuite
# 	@ ./testsuite

$(GEN)/%.o : %.cpp
	g++ $(CPPFLAGS) -c -o $@ $^

$(BIN)/buildRivers: $(GEN)/hydrologyFunctions.o $(GEN)/hydrologyParameters.o $(GEN)/buildRivers.o $(GEN)/hydrology.o $(GEN)/floatEndian.o $(GEN)/shore.o
	g++ $^ $(CPPFLAGS) -lopencv_imgproc -lgeos_c -o $@

$(BIN)/terrainPrimitives: $(GEN)/ts.o $(GEN)/terrainElevation.o $(GEN)/terrainHoneycomb.o $(GEN)/hydrology.o $(GEN)/floatEndian.o $(GEN)/terrainPrimitives.o $(GEN)/processTerrainPrimitives.o $(GEN)/shore.o
	g++ $^ $(CPPFLAGS) -lopencv_imgproc -lgeos_c -o $@

clean :
	@ rm -f $(BIN)/buildRivers $(BIN)/terrainPrimitives $(GEN)/*
#	@ cd $(GTEST_DIR)/scripts/test/ && make clean