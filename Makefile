# simplest possible to start.  dangerous since all header deps done manually. ne
rootPath = ./
include ${rootPath}/include.mk

all : vg2sg

clean : 
	rm -f  vg2sg
	cd sgExport && make clean
	cd tests && make clean

unitTests : vg2sg
	cd tests && make

vg2sg.o : vg2sg.cpp ${basicLibsDependencies}
	${cpp} ${cppflags} -I . vg2sg.cpp -c

${sgExportPath}/sgExport.a :
	cd ${sgExportPath} && make

vg2sg :  vg2sg.o ${basicLibsDependencies}
	${cpp} ${cppflags} ${basicLibs} vg2sg.o -o vg2sg 


test : unitTests
	pushd .  && cd ${sgExportPath} && make test && popd && tests/unitTests

