#Q Maintenance guide

####Update:Mauricio Esguerra March 26, 2014*
######Original by: John Marelius, August 29, 2000*


This document describes how to maintain the Q programs at the Åqvist
group.  Recently the code, documentation, and some scripts have been
added to version control using github, these documents still need some
updating. 


##File locations

The current source code for Q and documentation has been pushed to a
github private organization located at https://github.com/qusers


| What                                 | Where                                        |
|:------------------------------------ |:---------------------------------------------| 
| source code, makefile, history files | qsource/src /  history/                      |
| manual                               | qsource/documentation/manual/qman5.pdf       |
| license agreement                    | qsource/documentation/license.pdf            |
| this document                        | qsource/documentation/qmaintenanceguide.docx |
| force field files                    | qsource/ff/                                  |
| scripts                              | scripts/                                     |
| tests                                | qsource/tests/                               |
| Q web site                           | pending                                      |


##Updating Q

You need to be a member of the owners team at the github repository. 

###Updating source code on the web server

NOTE CHANGE THIS TO GIT WAY

~~After  modifying any  of  the master  source  files in  g:\\src\\q4:
Update  the  version number  and  last-modified-date  of  the file  in
question, e.g. md.f90, and in the main program, e.g. qdyn.f90. Run the
script g:\\src\\tarQ4.csh to update the files on the web server.~~

###Updating force field files

NOTE CHANGE THIS TO GIT WAY

~~After  modifying  the  master  files  in  g:\\FF4,  run  the  script
g:\\FF4\\ff2web.csh  to  update the  files  on  the  web server  (both
separated files and compressed archive files).~~

###Updating the manual

NOTE CHANGE THIS TO GIT WAY

~~After making  changes to the manual,  update the PDF  version on the
web server as follows: Log  on to Gem (outside Erling’s office).  Open
the manual  (g:\\doc\\Qman4.doc) Print it  to the print  driver called
‘Acrobat Distiller’. Check  the options to see that  Letter size paper
is used  (so that our US  users can print without  problem).  When the
Adobe Acrobat  window appears, save the  PDF file as  Qman4.pdf in the
wwwroot\\Q directory on the web server.~~



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
of programs except for the parallel version of qdyn5.

To compile using the Intel fortran compiler, and a debuggable version of the code
you would give the following options to the make command:

```bash
make debug COMP=ifort
```

###Windows

-   Not supported for now, we need a developer here.

###Linux - Generic(gfortran)
By default if a second COMP=[compilername] option is not given to the make command
the make program will use as defaults the gfortran compiler
```bash
make all
```

###Mac OSX (mavericks)
The compilation of the Q set of programs using the gfortran compiler needs a few different
options to those used in say, linux, so in order to compile in the most recent
version of the Mac OSX, you need to use:
```bash
make all COMP=osx
```

###Linux - CentOS 6.5 (triolith) Intel

```bash
module load intel/14.0.2
git clone https://github.com/qusers/qsource.git
cd qsource/src
cp makefile.ifort makefile
cp qdyn.F90_ifort_signals qdyn.f90
make all COMP=ifort
module load impi/4.1.3.048
make mpi COMP=ifort
```

###Linux - CentOS 6.5 (glenn) AMD Opteron 6220
There is some issue with finding the path to linker libraries, so, I've used
the compilation flags used in the mac in this case.
Needs benchmarking against ifortran and pgi.

```bash
git clone https://github.com/qusers/qsource.git
cd qsource/src
cp makefile.new makefile
module load  gcc/4.8/4.8.1
make all COMP=osx
module load openmpi/1.5.4
make mpi COMP=gcc
```                   


###Linux - CentOS 6.3 (csb) Intel

To compile at csb follow these steps:
```bash
source /home/apps/intel/composer_xe_2013.5.192/bin/ifortvars.sh intel64
export OMPI_FC=ifort
git clone https://github.com/qusers/qsource.git
cd qsource/src
cp makefile.new makefile
cp qdyn.F90_ifort_signals qdyn.f90
make all COMP=ifort
make mpi COMP=ifort
```

This will create the parallel executable qdyn5p

Once compilations is succesful you run the simple test in the /test folder
with  the run_test_mpi.sh script.  You  would only have to  create another script
for the slurm queue where you  make sure to load the openmpi libraries
with:
```bash
module load openmpi-x86_64
```

###Linux - Ubuntu 12.04 (abisko) AMD
To compile at abisko follow these steps:
```bash
module load openmpi/gcc/1.6.5
git clone https://github.com/qusers/qsource.git
cd qsource/src
cp makefile.new makefile
make all COMP=gcc

module load openmpi/gcc/1.6.5
make mpi COMP=gcc
```

###Linux - Scientific Linux 6.5 (tintin) AMD

To compile at tintin follow these steps:
```bash

unset SSH_ASKPASS

###USE YOUR USERNAME INSTEAD OF esguerra
git clone https://esguerra@github.com/qusers/qsource.git 

###It will ask for your password and download it showing something like:
remote: Counting objects: 934, done.
remote: Compressing objects: 100% (449/449), done.
remote: Total 934 (delta 476), reused 934 (delta 476)
Receiving objects: 100% (934/934), 16.79 MiB | 3.59 MiB/s, done.
Resolving deltas: 100% (476/476), done.

###If that works for you, hopefully, then go with:
cd qsource/src/
cp makefile.new makefile
module load gcc/4.8.2
make all COMP=gcc

###And to compile the parallel version
module load openmpi/1.4
make mpi COMP=gcc
```

##Updating the Q web site

The    Q    website    page    is    now   located    at    the    url
http://xray.bmc.uu.se/~aqwww/q, but  plans are underway  to migrate it
to a new address and keep it under version control.


##The E-mail list

The       home      page      of       the      list       is      at:
https://groups.google.com/d/forum/qmoldyn.     If   you're    in   the
google-group owners team  then you can see the  member list etc.  Just
log in  to your google account,  go to the  group web page and  on the
upper right side  of the page click on Manage. This  will show you all
the users  on the group (which  is configured as a  mailing list), and
also a lot of configuration information for the list.

To send a message to the list, address it to qmoldyn@googlegroups.com


