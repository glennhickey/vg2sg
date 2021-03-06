rootPath = ./
include ${rootPath}/include.mk

all : vg2sg

clean : 
	rm -f  vg2sg vglight.o pathspanner.o pathmapper.o vgsgsql.o vg2sg.o
	cd sgExport && make clean
	cd tests && make clean
	rm -f vg.pb.h vg.pb.cc vg.pb.o

unitTests : vg2sg
	cd tests && make

vg2sg.o : vg2sg.cpp vglight.h vg.pb.h ${basicLibsDependencies}
	${cpp} ${cppflags} -I . vg2sg.cpp -c

${sgExportPath}/sgExport.a : ${sgExportPath}/*.cpp ${sgExportPath}/*.h
	cd ${sgExportPath} && make

${protobufPath}/libprotobuf.a: ${protobufPath}/src/google/protobuf/*cc  ${protobufPath}/src/google/protobuf/*h
	cd ${protobufPath} && mkdir -p build && ./autogen.sh && ./configure --prefix=`pwd`/build/ && make && make install
	cp ${protobufPath}/build/lib/libprotobuf.a ${protobufPath}/

vg.pb.cc: vg.pb.h
vg.pb.h: vg.proto ${protobufPath}/libprotobuf.a
	${protobufPath}/build/bin/protoc vg.proto --cpp_out=.

vg.pb.o: vg.pb.h vg.pb.cc
	${cpp} ${cppflags} -I . vg.pb.cc -c 

vglight.o: vglight.cpp vglight.h vg.pb.h
	${cpp} ${cppflags} -I. vglight.cpp -c

pathmapper.o: pathmapper.cpp pathmapper.h pathspanner.h vglight.h vg.pb.h ${sgExportPath}/*.h
	${cpp} ${cppflags} -I. pathmapper.cpp -c

pathspanner.o: pathspanner.cpp pathspanner.h vglight.h vg.pb.h ${sgExportPath}/*.h
	${cpp} ${cppflags} -I. pathspanner.cpp -c

vgsgsql.o: vgsgsql.cpp vgsgsql.h pathmapper.h ${sgExportPath}/*.h
	${cpp} ${cppflags} -I. vgsgsql.cpp -c

vg2sg :  vg2sg.o vg.pb.o vglight.o pathspanner.o pathmapper.o vgsgsql.o ${basicLibsDependencies}
	${cpp} ${cppflags}  vg2sg.o vg.pb.o vglight.o pathspanner.o pathmapper.o vgsgsql.o  ${basicLibs} -o vg2sg 

test : unitTests
	cd ${sgExportPath} && make test && cd .. && tests/unitTests

