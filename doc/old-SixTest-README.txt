The SixDesk SixTrack test suite.  Eric.   16th May, 2014.

Just a few brief instructions for testing a new sixtrack version
interactively or on LSF batch at CERN. The procedures have been tested on
lxplus, lxbatch, MacOS and on Windows XP under cygwin.

NOTA BENE: SixTrack authors should submit a test case 
(normally fort.2 fort.3 fort.8 and fort.16) of a few 
thousand turns along with  the fort.6 and fort.10 output.

TESTS
The following test cases may be used (amongst others):
bnl exact crabamp javier thick4 lost tilt dipedge bb bb_ntwin1 notilt lostnotilt
lostevery s316 eric frs frs60 bbe51 bbe52 bbe571ib0 distance prob3 prob1 fma
elensidealthck4d elensidealthck6d elensidealthin4d elensidealthin6d.
To each directory there is a corresponding directory with suffix "in"
e.g.crabampin for crabamp containing the input files fort.2, 3, 8, and 16,
a fort.zip for BOINC, and fort.6.canonical and fort.10.canonical and 
fort.90.canonical.  for checking and comparing results. 
bnl[in] contains .dat input. There is no fort.10 for "exact".

The test cases to be used should be defined in the environment
variable TESTS e.g. setenv TESTS "crabamp javier" or
export TESTS="crabamp javier". The first few are rather short and can be
used for initial testing, the latter tests run for up to 100,000 turns.
(Other tests such as prob3 and prob1 run for up to one million turns.

There are three basic types of test: A) production with no Checkpoint/
Restart, B) production with one restart, and C) production with random kill
and restart. For tests B) and C) a Checkpoint/Restart version of SixTrack is
required.

EXECS
For each test type the variable EXECS should be defined as a list of the FULL paths 
to the SixTrack versions to be tested e.g.
setenv EXECS "/afs/cern.ch/user/m/mcintosh/scratch1/SixTracko0/SixTrack_4464_crlibm_bnl_ifort_O1/SixTrack_4464_crlibm_bnl_ifort_O1 \
/afs/cern.ch/user/m/mcintosh/scratch1/SixTracko0/SixTrack_4464_crlibm_bnl_ifort_sse2_O1/SixTrack_4464_crlibm_bnl_ifort_sse2_O1  \
/afs/cern.ch/user/m/mcintosh/scratch1/SixTracko0/SixTrack_4464_crlibm_bnl_ifort_sse3_O1/SixTrack_4464_crlibm_bnl_ifort_sse3_O1"
will test the three versions defined or will submit one batch job for each version.
For illustrative purposes and examples we shall assume two executables:
$HOME/SixTrack_9999_crlibm_bnl_ifort_sse2_O1  (a production version)
$HOME/SixTrack_9999_crlibm_bnl_ifort_boinc_api_sse2_O1 (a checkpoint/restart version)
and two tests crabamp and thick4.

NOTA BENE: BIN the directory try/bin contains some scripts and executables.
The executables used for checking are not (yet) available for MacOS but
the scripts should continue to run and the results can be SCP'd to lxplus
for comparison. On Windows CYGWIN you need to do:
 cd bin
 make clean
 make
which build with gfortran (likewise for macOS).

So now we define:
setenv TESTS crabamp thick4
setenv EXECS $HOME/SixTrack_9999_crlibm_bnl_ifort_sse2_O1

RUNNING the TESTS
It is worth noting that the following run scripts can be run
using e.g. the following sequence:
 sh
 nohup run > blog_run 2>&1 &
which allows longer tests to complete even after a logout.
The output in blog_run can be checked later.

A) run_pro [ID]
will run the two tests using $HOME/SixTrack_9999_crlibm_bnl_ifort_sse2_O1
and produce the results 
 fort.10.run_SixTrack_9999_crlibm_bnl_ifort_sse2_O1_[`hostname`]|ID]
 fort.6.run_SixTrack_9999_crlibm_bnl_ifort_sse2_O1_[`hostname`]|ID]
 fort.90.run_SixTrack_9999_crlibm_bnl_ifort_sse2_O1_[`hostname`]|ID]
(and fort.110_.... with most compilers). The bnl tests case is special
for the bnlelens tests and produces the *.dat_... files as well.
See later for result checking but the command
 comp_10 9999 
will compare the fort.10s with a subscript matching 9999 in the directories
$TESTS to the canonical results.

B) To test with a checkpoint/restart from the last checkpoint after all turns:
 setenv EXECS $HOME/SixTrack_9999_crlibm_bnl_ifort_boinc_api_sse2_O1
 run [ID]

C) and finally to run a test or two with random kills
 run_kill [ID]. 
Note that the ouput from these tests is rather verbose
due to set -x commands which are necessary for diagnosing
any errors. It can always just be deleted if noy required.

BATCH
To run tests on LSF batch the corresponding commands
 batch_pro [queue]
 batch [queue]
 batch_kill [queue]
can be used with results being returned to the $TESTS directories.
Note that the environment variables are inherited by the jobs.

The jobs are submitted to the queue 2nw by default or to the "queue" specified.
The batch scripts use the jobs batch_pro.lsf, batch.lsf and batch_kill.lsf
which in turn call the run scripts. The run scripts use corresponding
try scripts try_pro, try, and try_kill.

RESULTS
By default the script check_10 "pattern" compares the fort.10s in $TESTS
which match "pattern" to the fort.10.canonical (likewise check_90).
Similarly the script check_results will compare fort.6 as well but this is not
generally useful when different compilers are being used. The current canonical
versions were produced with ifort. The script comp_10 is similar to check10
but in addition prints the magnitude of result differences.
Finally, the script check_extras will compare any files listed in
the text file extra_checks.txt in the test folder using
the tool (e.g. "diff") listed.

Note that fort.10's should be IDENTICAL (the comparisons ignore the SixTrack version
Word 52 of each line in fort.10 and the CPU time in Word 60).
The results from ifort opt 0 or 1, ia32, sse2, sse3, nagfor, gfortran, pgf90 and lf95 (opt 0)
have been checked and are identical (and fort.90 too of course).
