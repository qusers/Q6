#Q Maintenance guide

####Update:Mauricio Esguerra October 20, 2015*
######Original by: John Marelius, August 29, 2000*


This document describes  how to maintain the Q programs  at the Åqvist
lab.  
In  2014 the  code,  documentation,  and some  scripts  were added  to
version  control  using  github,   note that these  documents  need  constant
updating.


##File locations

The current source code for **Q** and documentation has been pushed to
a github private organization located at https://github.com/qusers


| What                                 | Where                                        |
|:------------------------------------ |:---------------------------------------------| 
| source code, makefile, history files | qsource/src / qsource/history/               |
| manual                               | qsource/documentation/manual/qman.pdf        |
| license agreement                    | qsource/documentation/license.pdf            |
| this document                        | qsource/documentation/qmaintenanceguide.md   |
| force field files                    | qsource/ff/                                  |
| scripts                              | scripts/                                     |
| tests                                | qsource/tests/                               |
| Q web site                           | http://qdyn.no-ip.org                        |


##Updating Q

You need to be a member of the owners team at the github repository.


###Updating source code on the web server

After modifying  any of  the master source  files with  .f90 extension
update  the  version   number  (QDYN_VERSION)  and  last-modified-date
(QDYN_DATE) of  the file  in question,  e.g.  md.f90, and  in the  main
program, e.g. qdyn.f90.


###Updating force field files

Force-field files are stored in  the ff folder. Github versioning will
take care  of changes  to force-field  parameters (.par)  and topology
(.lib) files.


###Updating the manual

The **Q** manual is kept at documentation/manuals/qman.tex. Since it's
a LaTeX plain  text file all changes are versioned  directly by git. A
current version  of the  pdf manual  should be  uploaded to  the **Q**
webserver when changes are made.


##Building executables

This sections describes how to compile the code using the makefile to
select which pakages to compile and which compiler.  
In short, a new makefile updated to use te GNU fortran compilers, and
the intel fortran compiler has been created.  

In general the use of the makefile is as follows:  

```bash
make
```

The make alone command will show you the options which have to be specified at
compilation depending on the compiler available in your machine and the set of
programs you want to compile, for example, the [all] option stands for compilation
of programs except for the parallel version of qdyn.  

To compile  using the  Intel fortran compiler,  and a  debuggable (via
gdb) version of  the code you would give the  following options to the
make command:  

```bash
make debug COMP=ifort
```


###Windows

Compiling in Windows XP has been tested to work for the serial version
of the code, not the parallel one yet.  

The successful compilation was done using  MinGW to be able to use the
GCC  (Gnu  Compiler Collection),  therefore  you  will  have to  first
install MinGW, then make sure that  the location of it is in your path
(C:\MinGW\bin) and then  you can compile in the  command prompt in the
same way it's done in unix flavors, that is:  

```bat
make all COMP=gcc
```

It  is quite  possible that  the  compilation will  go trough  without
issues with the intel compiler, but this has not been tested yet.


###Linux - Generic(gfortran)

By default if a second  COMP=[compilername] option is not given to the
make  command the  make  program  will use  as  defaults the  gfortran
compiler options:  

```bash
make all
```


###Mac OSX (mavericks)

The  compilation of  the  **Q**  set of  programs  using the  gfortran
compiler needs a few different options to those used in say, linux, so
in order  to compile in  the most recent version  of the Mac  OSX, you
need to use:

```bash make all COMP=osx ```


###Linux - CentOS 6.5 (triolith) Intel

```bash
module load intel/14.0.2
git clone https://github.com/qusers/qsource.git
cd qsource/src
make all COMP=ifort
module load impi/4.1.3.048
make mpi COMP=ifort
```


###Linux - CentOS 6.6 (glenn) AMD Opteron 6220  

For hardware info on the cluster go to:

    http://www.c3se.chalmers.se/index.php/Hardware_Glenn

```bash
unset SSH_ASKPASS
git clone https://yourusername@github.com/qusers/qsource.git
cd qsource/src
module load intel-compilers/15.0/090
make all COMP=ifort
module load intel-mpi/5.0.1.035
make mpi COMP=ifort
```

After compilation make sure to run the benchmarks changing the options
to  adapt to  the  slurm manager  at  Glenn and  also  make sure  path
declarations point to your **Q** binaries.


###Linux - CentOS 6.3 (csb) Intel

To compile at csb follow these steps:  

```bash
source /home/apps/intel/composer_xe_2013.5.192/bin/ifortvars.sh intel64
export OMPI_FC=ifort
git clone https://github.com/qusers/qsource.git
cd qsource/src
make all COMP=ifort
make mpi COMP=ifort
```

This will create the parallel executable qdynp  

Once compilations  is succesful you run  the simple test in  the /test
folder with the run_test_mpi.sh script.  You would only have to create
another script for the slurm queue where you make sure to load the
openmpi libraries with:  

```bash
module load openmpi-x86_64
```


###Linux - Ubuntu 12.04 (abisko) AMD

To compile at abisko follow these steps:  

```bash
module load openmpi/gcc/1.6.5
git clone https://github.com/qusers/qsource.git
cd qsource/src
make all COMP=gcc
module load openmpi/gcc/1.6.5
make mpi COMP=gcc
```


###Linux - Scientific Linux 6.7 (tintin) AMD

At  tintin both  gcc and  intel fotran  compilers are  available.  The
architecture of  the computer nodes is  of the AMD family.   Each node
has two Opteron 6220 cpu's having 8 nodes each.  

To compile at tintin using gcc use this recipe:  

```bash

unset SSH_ASKPASS

###USE YOUR USERNAME INSTEAD OF esguerra
git clone https://esguerra@github.com/qusers/qsource.git 

###It will ask for your password and download the full qsource repository 
###containing code and documentation showing something like:
remote: Counting objects: 934, done.
remote: Compressing objects: 100% (449/449), done.
remote: Total 934 (delta 476), reused 934 (delta 476)
Receiving objects: 100% (934/934), 16.79 MiB | 3.59 MiB/s, done.
Resolving deltas: 100% (476/476), done.

###If that works for you, hopefully, then go with:
cd qsource/src/
module load gcc/4.8.2
make all COMP=gcc

###And to compile the parallel version
module load openmpi/1.4
make mpi COMP=gcc
```

For compilation using the ifortran compiler:  

```bash
module load intel openmpi
make all COMP=ifort
make mpi COMP=ifort
```


##Updating the Q web site

The    Q    website    page    is    now   located    at    the    url
http://xray.bmc.uu.se/~aqwww/q, but  plans are underway  to migrate it
to a new address and keep it under version control. For now the page 
can also be found at the tempory address:  

    http://qdyn.no-ip.org/  

Where      versioning      is       being      used      from      the
https://github.com/qusers/qwebsite repository.


##The E-mail list

The       home      page      of       the      list       is      at:
https://groups.google.com/d/forum/qmoldyn.     If   you're    in   the
google-group owners team  then you can see the  member list etc.  Just
log-in  to your google account,  go to the  group web page and  on the
upper right side  of the page click on Manage. This  will show you all
the users  on the group (which  is configured as a  mailing list), and
also a lot of configuration information for the list.

To send a message to the list, address it to qmoldyn@googlegroups.com


