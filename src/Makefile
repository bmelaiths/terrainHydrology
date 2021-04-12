CPPFLAGS += -g -Wall

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

test : kdtreetest
	@ ./$<

kdtree.o : kdtree.cpp kdtree.hpp

buildRivers.o : buildRivers.cpp kdtree.hpp

buildRivers: buildRivers.o kdtree.o test
	g++ buildRivers.o kdtree.o -lopencv_core -lopencv_imgproc -lopencv_highgui -o $@

clean :
	@ rm -vf buildRivers kdtreetest *.gch *.o
	@ cd $(GTEST_DIR)/scripts/test/ && make clean