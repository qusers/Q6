Q Maintenance guide

*John Marelius, August 29, 2000*

\

This document describes how to maintain the Q programs on the Åqvist
group computers. 

File locations

The locations of the Q files are listed in the table below, where it is
assumed that device g: is mapped to \\\\skalleper\\bin and device w: to
\\\\skalleper\\inetpub\$.

**whatwhere**source code, makefile, Windows project files, history
filesg:\\src\\q4manualg:\\doc\\Qman4.doclicence
agreementg:\\doc\\Q\_license.docthis documentg:\\doc\\Q Maintenance
guide.docforce field filesg:\\FF4samplesg:\\qsamplesQ web site -
editable copyw:\\wwwdev\\QQ web site - files read by web
serverw:\\wwwroot\\QUpdating Q

You need to be a member of the NT user group WWW Admins to have
permission to do some of these steps. Note that the procedures below
will update files in the wwwroot directory of the web server directly
rather than updating the copies in the wwwdev directory. Take care not
to overwrite files in wwwroot with older versions from wwwdev!

Updating source code on the web server

After modifying any of the master source files in g:\\src\\q4:

Update the version number and last-modified-date of the file in
question, e.g. md.f90, and in the main program, e.g. qdyn.f90.

Run the script g:\\src\\tarQ4.csh to update the files on the web server.

Updating force field files

After modifying the master files in g:\\FF4, run the script
g:\\FF4\\ff2web.csh to update the files on the web server (both
separated files and compressed archive files).

Updating the manual

After making changes to the manual, update the PDB version on the web
server as follows:

Log on to Gem (outside Erling’s office).

Open the manual (g:\\doc\\Qman4.doc)

Print it to the print driver called ‘Acrobat Distiller’. Check the
options to see that Letter size paper is used (so that our US users can
print without problem).

When the Adobe Acrobat window appears, save the PDF file as Qman4.pdf in
the wwwroot\\Q directory on the web server.

Building executables

This section describes how to transfer source code, compile and package
the executables.

Windows NT/Alpha

Log on to one of the NT Alphas.

In Visual Fortran, use the batch build command to build all the
‘Release’ executables (optimisation on), i.e. the configurations
Qcalc4‑AlphaRel, Qdyn4‑ AlphaRel, Qfep4‑ AlphaRel, Qprep4‑ AlphaRel and
Qdyn4‑AlphaDum.

Run the script g:\\src\\ ZipQ4bin\_x86.csh to make a zip file on the web
server with the executables.

\

Windows NT/x86

Log on to Gem (outside Erling’s office).

In Visual Fortran, use the batch build command to build all the
‘Release’ executables (optimisation on), i.e. the configurations
Qcalc4‑IntelRel, Qdyn4‑IntelRel, Qfep4‑IntelRel, Qprep4‑IntelRel and
Qdyn4‑IntelDum.

If required, build also the evaluation executables (optimisation off),
i.e. the configurations Qcalc4‑IntelEval, Qdyn4-IntelEval,
Qfep4-IntelEval and Qprep4‑IntelEval.

Log on to one of the NT Alphas.

Run the script g:\\src\\ ZipQ3bin\_x86.csh to make a zip file on the web
server with the executables.

If required, also run the script g:\\src\\ zipq4eval.csh to make a zip
file with the evaluation kit.

Compaq UNIX

Transfer the whole source code archive (update it first!) by ftp to the
directory /aq/sudd/john/src/q4 on dqs1.

Log on to dqs1.

Unpack the source code archive:\
 gunzip -c Q4\_src.tar.gz | tar -xvf -

Alternatively, replace steps 1 and 2 by transferring individual source
files.

Remove old object files:\
 make clean

Build the executables for use on the dqs cluster (optimised for the
Alpha EV6 or 21264 processor):\
 make alpha-osf1-ev6

To update the locally used executables, copy them to the shared
executable directory:\
 cp Q\*4 \$BIN/

Package the executables:\
 package\_ev6.csh

Remove object files\
 make clean

Build the executables for use on the other UNIX Alpha workstations
(optimised for the Alpha EV5 or 21164 processor):\
 make alpha-osf1-ev5

To update the locally used executables, copy them to the shared
executable directory:\
 cp Q\*4 /aq/sudd/john/bin/ OSF1\_V4.0\_EV5 /

Package the executables:\
 package\_ev5.csh

Retrieve the archive files Q4\_DigitalUNIX4\_Alpha21264.tar.gz and
Q4\_DigitalUNIX4\_Alpha21164.tar.gz by ftp to the directory
wwwroot\\Q\\download\\bin on the web server.

SGI Irix

Transfer the whole source code archive (update it first!) by ftp to a
directory named q4 in your home directory on genome.ibg.uu.se.

Log on to genome.ibg.uu.se with SSH.

Unpack the source code archive:\
 gunzip -c Q4\_src.tar.gz | tar -xvf -

Alternatively, replace steps 1 and 2 by transferring individual source
files.

Remove old object files:\
 make clean

Build the executables for use on MIPS R10000 processors:\
 make irix6.4

Package the executables:\
 ./package\_R10k.csh (copy this script from ˜john/src/q4)

Remove object files\
 make clean

Build the executables for use on MIPS R4000 processors:\
 make irix6.2

Package the executables:\
 ./package\_R4k.csh (copy this script from ˜john/src/q4)

Retrieve the archive files Q4\_IRIX6.4\_R10k.tar.gz and
Q4\_IRIX6.2\_R4k.tar.gz by ftp to the directory
wwwroot\\Q\\download\\bin on the web server.

Linux

Ask Peter Vagedes \<vagedes@chemie.fu-berlin.de\> to build the
executables.

Updating the Q web site

The Åqvist group web site http://aqvist.bmc.uu.se and the Q web site
http://aqvist.bmc.uu.se/Q can both be edited using VisualPage. The
VisualPage project file is w:\\aqvist\_dev.vpp. The HTML and other files
in this project are all in the w:\\wwwdev directory. Edit the files
there and then copy the modified files, or all files, to w:\\wwwroot
which is the actual web server directory.

The E-mail list

The home page of the list is http://Qusers.listbot.com. To see the
member list etc, log in as list owner using the E-mail address
John.Marelius@molbio.uu.se and the password agent007

To send a message to the list, address it to Qusers@listbot.com
