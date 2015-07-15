binPath=${rootPath}
libPath=${rootPath}

#Modify this variable to set the location of sonLib
vgRootPath=${rootPath}/../vg
vgLibPath=${vgRootPath}

sgExportPath=${rootPath}/sgExport

include ${sgExportPath}/include.mk

cflags += -I ${vgLibPath} -I ${sgExportPath}
cppflags += -I ${vgLibPath} -I ${sgExportPath}
basicLibs = ${vgLibPath}/libvg.a ${sgExportPath}/sgExport.a 
basicLibsDependencies = ${basicLibs}


