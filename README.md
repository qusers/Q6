Q6
=======
Version control of the molecular dynamics code called Q, version 6, via github.

Q is a set of Molecular Dynamics (MD) tools tailored to the following specific kinds of free energy calculations:

1. Free Energy Perturbation (FEP)
2. Empirical Valence Bond (EVB)
3. Linear Interaction Energies (LIE)
4. Quantum Classical Path (QCP)

© Copyright 2017 Johan Åqvist, John Marelius, Shina Caroline Lynn Kamerlin and Paul Bauer.

Developed by:	The Department of Cell and Molecular Biology
	        Uppsala University, Uppsala, Sweden
	        http://xray.bmc.uu.se/~aqwww/q/ 
		paul.bauer.q@gmail.com, qmoldyn@googlegroups.com 

## Installation
The current makefiles make it relatively easy to compile the code in a linux or mac OSX environment.

You have to make sure first that gfortran is installed and if you want to compile the MPI version
you also have to make sure that openMPI is installed.

If you're using gfortran make sure that you have version 4.8 or later. This applies both to compilation in mac and linux.
To make sure that you have version 4.8 or later use:

```bash
gfortran --version
```

Right before issuing the "make" command.

To install in a linux environment you have to move to the src/ folder where the source code and the makefile are located at. To get information on how to use the makefile just type **make** in your terminal:

```bash
unset SSH_ASKPASS
git clone https://www.github.com/qusers/Q6.git
cd Q6/src
make
```
Following the instructions for compilation using the gfortran compiler and the non-mpi version of the code would then take the form:

```bash
make all COMP=gcc
```

Or for MPI compilation (after loading a proper MPI compiler library):
```bash
make mpi COMP=gcc
```

To compile in Mac OSX 10.9.2 you can use native gfortran binaries which you can download from  (https://gcc.gnu.org/wiki/GFortranBinaries#MacOS) or you can also compile using the GCC (Gnu Compiler Collection) distributions available to fink, macports or homebrew. In order to compile in a MAC you should call the fink environment first. Usually with:

```bash
source /sw/bin/init.sh
```

alternatively, you can install [homebrew](http://brew.sh/) and use (confirmed for Mac OSX 10.10):
```bash
brew install gcc
```


Then you proceed as in linux
```bash
git clone https://www.github.com/qusers/Q6.git
cd Q6/src
make all COMP=osx
```

This will take care of redirecting the binaries and object files to standard bin and obj folders.

After this you have to add the program to your system path by modifying your shell initiation script, that is, if your shell is bash, you can add the following lines to your .bashrc file using a text editor:

```bash
SOFT=/Users/johndoe/software
export QDIR=$SOFT/Q6
export PATH=$QDIR/bin:$PATH  
```
Where $SOFT will be the place where your software folder is located at, e.g. /Users/johndoe/software

Once the Q6 binaries are declared in your path you should be able to call all Q6 binaries from your terminal.
To test that the path to your compiled Q binaries has been correctly assigned you can issue the following commands in the terminal:

```bash
source .bashrc
echo $path | grep Q6

/Users/johndoe/software/Q6

Qprep6

Build and version information

Build number 6.0.X
Build date   20XXXXXX
Built:       
      by     johndoe
      on     localhost
      git id 961bca2e036ca56218d1cca7a134472fd7e3cc86
      with   GNU Fortran Debian 6.4.0-1 6.4.0 20XXXXXX
Qprep version 6.0.X , git id=961bca2e036ca56218d1cca7a134472fd7e3cc86 initialising
Current date 20XX-XX-XX and time XX:XX:XX

Welcome in Qprep modification date 20XXXXXX

 Q6, Copyright © 2017 Johan Åqvist, John Marelius, Shina Caroline Lynn Kamerlin and Paul Bauer
 Q6 comes with ABSOLUTELY NO WARRANTY.  This is free software, and you are welcome
 to redistribute it under certain conditions. For details, add the --help flag.

Qprep> 
```

