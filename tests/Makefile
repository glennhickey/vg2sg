rootPath = ..
include ${rootPath}/include.mk

all : unitTests

clean :
	rm -f *.o unitTests

unitTests : CuTest.o unitTests.o pathmapperTests.o 

CuTest.o : ${sgExportPath}/tests/CuTest.h ${sgExportPath}/tests/CuTest.c
	${cxx} ${cflags} -I ${sgExportPath}/tests -c ${sgExportPath}/tests/CuTest.c	

unitTests.o : unitTests.h unitTests.cpp
	${cpp} ${cppflags} -I ${sgExportPath}/tests -I${rootPath}/ unitTests.cpp -c

pathmapperTests.o : pathmapperTests.cpp ${rootPath}/*.h
	${cpp} ${cppflags} -I ${sgExportPath}/tests -I${rootPath}/ pathmapperTests.cpp -c

unitTests : CuTest.o unitTests.o pathmapperTests.o ${sgExportPath}/sgExport.a
	${cpp} -I ${sgExportPath}/tests -I${rootPath}/ ${cppflags} ${sgExportPath}/sgExport.a CuTest.o unitTests.o pathmapperTests.o  -o unitTests

