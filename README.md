VOFI Library
======

The  VOFI library initializes the volume fraction scalar field 
in a computational mesh with cubic cells given an analytic expression 
f(x,y,z) for the interface. The implicit function f is specified
by the user, the interface is the zero level set of the function,
f(x,y,z) = 0, and the reference phase is located where f(x,y,z) < 0. 
Each routine in the directory 'src' contains a brief description
of what it does and of the I/O variables. 

The algorithm implemented in the library has been described in the following papers:

[1] S Bnà, S Manservisi, R Scardovelli, P Yecko, S Zaleski, "Numerical integration of implicit functions for the initialization of the VOF function", Computers & Fluids 113, 42-52, [doi:10.1016/j.compfluid.2014.04.010] 

[2] S Bnà, S Manservisi, R Scardovelli, P Yecko, S Zaleski, "VOFI -- A library to initialize the volume fraction scalar field", Computer Physics Communications, Computer Physics Communications, 2015, [doi:10.1016/j.cpc.2015.10.026]

## Vofi-specific configuration options

#### Build Process

* Extract vofi-<package_version>.tar.gz archive into any temporary directory

* Choose a local (e.g. /home/username) or system (e.g. /usr/local)
    destination for '/installing_directory' which is write-able

* The default c-compiler flags are: 
    "-O2 -fomit-frame-pointer -ffast-math -Wall". 
    The user can change the default value by setting the environment 
    variable CFLAGS 
    (e.g. "export CFLAGS=-g" to compile the library in debug mode)
    
* ./configure --prefix=/installing_directory
    (by default, a shared library is built on platforms that support it;
     the user can specify modified forms of the configure flags 
     "--enable-shared" and "--enable-static" to choose whether shared 
     and/or static libraries are built)
     
* make all
    (to build the library and to compile tests) 

    if separately:

    cd src
    make all 
    (to build only the library) 

    cd demo_src
    make all
    (to compile all C, Fortran and CPP tests)

    cd demo_src/C
    make all
    (to compile only C tests)
    
* make check
    (in demo_src: run all C, Fortran and CPP tests in the temporary 
     directory; in demo_src/C: run only C tests.
     A tests uite outcome for C, CPP and Fortran tests is printed both on 
     terminal and in the file test-suite.log, containing only the summary 
     of the exit status of each executable.
     A log file is also generated for each test in the C, CPP and Fortran 
     folders where the user can check the integration result in his/her 
     platform and in our computer, for example the file cap1_c.log) 
     
* make install
    (if the installing_directory is /usr/local, then libvofi.a, 
     or libvofi.so, will be located in /usr/local/lib64, vofi.h in 
     /usr/local/include and the demo executables in /usr/local/bin)
     
     
* make installcheck
    (run all C, Fortran and CPP tests with the installed version of the 
     library. The output of each test is printed on terminal together 
     with the results obtained in our computer )

* a few routines can run with random input parameters with a seed 
    based on the CPU clock, in demo_src/C these are the following 
    five executables: 
    ellipse_c, gaussian_c, rectangle_c, sine_line_c, sine_surf_c. 

    To run the ellipse test with static parameters type 
    ./ellipse_c

    while with random parameters type one of the two
    ./ellipse_c -r 
    ./ellipse_c --randominput

    The same approach applies also to the corresponding executables in 
    demo_src/CPP and demo_src/Fortran
    
    
    
## Requirements

The vofi library has no particular requirement other than
reasonably modern C, Fortran and C++ compilers. 


## Brief description of the files and directories in the distribution

#### Files in the top directory

* configure.ac: file to generate the configure script by autoconf
* makefile.am:  file to generate makefile.in by automake
* README:       this file

#### The other files

    aclocal.m4    ChangeLog     config.guess   config.sub 
    configure     depcomp       INSTALL        install-sh 
    ltmain.sh     makefile.in   missing        NEWS
    test-driver   AUTHORS       COPYING

have been generated automatically by the GNU autotools. 


#### Subirectories of the top directory

    demo_src   include   m4   src
    
#### Subdirectory demo_src

It contains the source files for the tests. At the first two levels 
there is also the file makefile.am to generate makefile.in by automake    

There are three subdirectories: C, CPP, Fortran

Each of them contains two more subdirectories: 2D, 3D

The directory '2D' contains the two-dimensional tests in the 
subdirectories:

Ellipse   Gaussian   Rectangle   Sine_line 

The directory '3D' contains the three-dimensional tests in the 
subdirectories:

Cap1   Cap2   Cap3   Sine_surface   Sphere 

All the directories with a test contain a main program and a file 
with the analytical expression for the line or the surface, for tests 
in C/C++ there is also an include file

The user can change the input parameters at will, but should also 
change the value of the area/volume accordingly, and recompile the tests 

As an example, we consider the test for an ellipse in C, with the 
following relative path from the temporary directory:

    ./vofi-1.0/demo_src/C/2D/Ellipse

with the three files: ellipse.c  ellipse.h  main_ellipse.c

* ellipse.c: it contains the implicit function 'impl_func' with the 
           analytical expression of an ellipse and the function 
           'check_area' to compute the ellipse area and to compare this 
           value with that obtained in our computer

* ellipse.h: input parameters for the grid resolution, the computational 
           box, position of the ellipse center, its semi-axes and their 
           angle with respect to the Cartesian axes

* main_ellipse.c: it defines the grid and calls the functions to initialize 
                the volume fraction and to check the area 

#### Subdirectory include:

It contains three include files:

    vofi_GL.h   vofi_stddecl.h   vofi.h

* vofi_GL.h: it contains nodes and weights for the Gauss-Legendre's integration

* vofi_stddecl.h: it contains several declarations used by the library routines

* vofi.h: it contains the functions prototype and should be included by
        the user when calling the VOFI routines  
        

#### Subdirectory m4:

It contains the Libtool macro files:

    libtool.m4   lt~obsolete.m4   ltoptions.m4   ltsugar.m4   ltversion.m4

These macro processor files are used toghther with aclocal.m4 by the GNU 
autotools (aclocal, libtoolize and autoconf) to generate the configure 
script (configure). During configuration, the libtool script is generated 
by config.status and ltmain.sh files. The libtool is then invoked as needed 
at run time by the make command as a wrapper around compilers, linkers, 
install and cleanup programs.      


#### Subdirectory src:

Besides the two files makefile.am and makefile.in, it contains
the ten source files of the library:

    checkconsistency.c   getcc.c       getdirs.c   getfh.c
    getintersections.c   getlimits.c   getmin.c    getzero.c
    integrate.c          interface.c
        

* checkconsistency.c: it contains two functions to check the consistency
                    with a minimum on a cell side and on a cell face;
                    this greatly reduces the number of calls to the
                    functions that compute a minimum


* getcc.c: driver to compute the integration limits and the volume fraction 
         in two and three dimensions 


* getdirs.c: it checks if the cells is either full or empty, if not
           it determines the main, second and third coordinate directions


* getfh.c: it computes the characteristic function value fh


* getintersections.c: it contains two functions to compute the interface
                    intersection(s) with a cell side and inside a face,
                    these are internal/external limits of integration


* getlimits.c: it subdivides the side along the secondary or tertiary
             direction to define rectangles or rectangular hexahedra 
             with or without the interface 


* getmin.c: it contains two functions to compute the function minimum 
          either in a given segment or in a cell face, the search is
          stopped if a sign change is detected 


* getzero.c: it computes the zero in a given segment 


* integrate.c: it contains two functions to compute the normalized cut 
             area/volume with a single/double Gauss-Legendre quadrature 


* interface.c: it contains two functions to call from Fortran the
             corresponding C functions
