Compilation notes:
 
- Compilation works on Unix-like systems (Linux, MacOSX, Cygwin)
- you need g++, make, cmake (command line tools) installed on your system
- if you are under Cygwin, you need to replace the file make.exe in C:\cygwin\bin by this one:
    http://www.cmake.org/files/cygwin/make.exe

To compile, either simply execute scriptcompile in a terminal from this folder, or go to cmake directory and type the following commands:

cmake .. -G"Unix Makefiles" -DCMAKE_BUILD_TYPE:STRING=Release
make

The binaries will be placed in the "bin" folder

Optional:
- by default the code is configured for 3-dimensional data. You can change this setting by adding
    "-DDIM=n" at the end of the cmake command (where n is the dimension you want)
- by default the compilation process checks if Cuda is installed on your system. If you want to
    disable Cuda, add "-DENABLECUDA=NO" at the end of the cmake command above. Compilation with Cuda on Cygwin
    does not work currently, so Cuda will be disabled on this platform, whatever setting you choose.


    



