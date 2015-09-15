#copied from sonLib...

# we do specific stuff for specific host for now.
HOSTNAME = $(shell hostname)
MACH = $(shell uname -m)
SYS =  $(shell uname -s)

cxx = gcc -std=c99
cpp = g++ 

# Compiler flags.
# DO NOT put static library -l options here. Those must be specified *after*
# linker input files. See <http://stackoverflow.com/a/8266512/402891>.

#Release compiler flags
cflags_opt = -O3 -g -Wall --pedantic -funroll-loops -DNDEBUG 
#-fopenmp
cppflags_opt = -O3 -g -Wall -funroll-loops -DNDEBUG -Wno-sign-compare

#Debug flags (slow)
cflags_dbg = -Wall -Werror --pedantic -g -fno-inline
cppflags_dbg = -Wall -g -O0 -fno-inline 

#Ultra Debug flags (really slow)
cflags_ultraDbg = -Wall -Werror --pedantic -g -fno-inline

#Profile flags
cflags_prof = -Wall -Werror --pedantic -pg -O3 -g

#for cpp code: don't use pedantic, or Werror
cppflags = ${cppflags_opt} 

#Flags to use
cflags = ${cflags_dbg} 


binPath=${rootPath}
libPath=${rootPath}

sgExportPath=${rootPath}/sgExport
protobufPath=${rootPath}/protobuf

cflags +=  -I ${sgExportPath}
cppflags +=  -I ${sgExportPath} -I ${protobufPath}/build/include
basicLibs = ${sgExportPath}/sgExport.a ${protobufPath}/libprotobuf.a -lz
basicLibsDependencies = ${sgExportPath}/sgExport.a ${protobufPath}/libprotobuf.a


