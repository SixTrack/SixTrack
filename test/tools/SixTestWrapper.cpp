/*
// This prevents building on SL6
#ifdef _GNU_SOURCE
#undef _GNU_SOURCE
#endif
*/

#include <cstdlib>
#include <iostream>
#include <fstream>

#include <string>
#include <sstream>
#include <vector>
#include <list>
#include <algorithm>

#include <unistd.h>
#include <signal.h>

#ifndef WIN32
//open()
#include <fcntl.h>
#include <sys/stat.h>
//waitpid()
#include <sys/types.h>
#include <sys/wait.h>
#endif

//to get errors from POSIX calls
#include <errno.h>
#include <stdio.h>

//memcmp
#include <string.h>

//pthread
#include <pthread.h>

#ifdef WIN32
#include <windows.h>
#endif

#ifdef LIBARCHIVE
#include <sys/stat.h>
#include <dirent.h>
#include "libArchive_wrapper.h"
#endif

#ifndef WIN32
void *pthread_kill_sixtrack(void*);
void *pthread_wait_sixtrack(void*);
#else
DWORD WINAPI winthread_kill_sixtrack(LPVOID);
DWORD WINAPI winthread_wait_sixtrack(LPVOID);
#endif

bool CopyFile(std::string InputFileName, std::string OutputFileName);
bool FileComparison(std::string f1, std::string f2);
size_t StripCR(std::string);

void CheckPrint(std::string,bool);
void PrettyDivider(std::string);
bool CheckFort10(char**);
bool CheckFort90(char**);
bool CheckSTF(char**);
bool PerformExtraChecks(bool&, char* convert_dump_bin, char* dump_bin_files);
std::vector<int> ParseKillTimes(char*);

void UnlinkCRFiles();

#ifndef WIN32
struct KillInfo {
    pid_t SixPID;       //The pid
    int   kTime;        //The time to kill for
    bool  RunStatus;    //Did the run finish whilst we were waiting to kill?
    bool  CRKILLSWITCH; //SixTrack was stopped by the internal CRKILLSWITCH
};
#else
struct KillInfo {
    HANDLE SixHANDLE;   //The handle
    DWORD kTime;        //The time to kill for
    BOOL RunStatus;     //Did the run finish whilst we were waiting to kill?
    BOOL CRKILLSWITCH;  //SixTrack was stopped by the internal CRKILLSWITCH
};
#endif

/**
* SixTrack testing wrapper
* This runs a given SixTrack binary (given as the first argument)
* It then performs checks on the output
* Currently we check fort.10 and fort.90
* Arguments:
* 1: Path to Sixtrack binary to run
* 2: Path to checkf10 binary to run
* 3: Path to read90 binary to run
* 4: bool to check fort.10
* 5: bool to check fort.90
* 6: bool for STF enabled
* 7: Number of files to expect in Sixout.zip (0 means no Sixout.zip)
* 8: CR enabled
* 9: CR kill time
* 10: Path to readDump3 binary to run
* 11: Name of files (comma-separated list) in extra_checks that is a dump format 3/8 (binary) or NONE
* 12: Minimum number of sucessfull checkpointed restarts to pass the test [int]
*
* For running the tools:
* On "Unix" we call fork() and exec()
* On windows we call spawn()
* https://docs.microsoft.com/en-us/cpp/c-runtime-library/spawn-wspawn-functions
* https://docs.microsoft.com/en-us/cpp/c-runtime-library/reference/spawnl-wspawnl
*
*/
int main(int argc, char* argv[])
{
    //First check we have the correct number of arguments
    if(argc != 13) {
        std::cout << argv[0]
                  << " called with the incorrect number of arguments, should be 12, "
                  << "but was called with " << argc - 1 << " arguments" << std::endl;
        return EXIT_FAILURE;
    }

    // Want our first argument to be the sixtrack binary name, assume this is run in the test folder.
    std::cout << std::endl;
    std::cout << "SixTestWrapper called with " << argc << " arguments:" << std::endl;
    for(int n=0; n < argc; n++) {
        printf("%4d : %s\n",n,argv[n]);
    }
    std::cout << std::endl;

    //Set defaults:
    //Set to true if we have a CR build
    bool CRon = false;

    //How long to wait in seconds before killing the CR run
    std::vector<int> KillTimes;
    int KillTime = 0;

    bool fort6 = false;
    bool fort10 = false;
    bool fort90 = false;
    bool STF = false;
    bool extrachecks = false; //Returned from PerformExtraChecks
    int  sixoutzip = 0;

    bool fort10fail = false;
    bool fort90fail = false;
    bool STFfail = false;
    bool ExtraChecksfail = false;
    bool sixoutzipfail = false;
    bool crRestartFail = false;

#ifdef LIBARCHIVE
    const char* const tmpdir = "sixoutzip_tmpdir";
    const char* const sixoutzip_fname = "Sixout.zip";
#endif

    if(atoi(argv[4]) != 0) {
        fort10 = true;
    }

    if(atoi(argv[5]) != 0) {
        fort90 = true;
    }

    if(atoi(argv[6]) != 0) {
        STF = true;
    }
    if(atoi(argv[7]) != 0) {
        sixoutzip = atoi(argv[7]);
    }

    if(atoi(argv[8]) != 0) {
        CRon = true;
        KillTimes = ParseKillTimes(argv[9]);
    }

    /**
     * First step is to handle all the running of sixtrack.
     * This include Checkpoint/resume builds where the run will be killed after a set period of time and then restarted.
     * In each case, to run sixtrack we call fork and let the child process do the exec()
     * The main thread then calls waitpid() for sixtrack to finish executing.
     * For CR we wait then send a kill signal after a specified number of seconds
     *
     * The following code is *really* ugly.
     *
     */
    if(CRon) {
        std::cout << "  Starting CR run loop - will clear out any files from a previous run" << std::endl;
        UnlinkCRFiles();

        size_t KillCount=0;
        while(true) {

            if (KillCount >= KillTimes.size()){
                std::cout << "  Reached maximum KillCount = " << KillCount << ", exiting now." << std::endl;
                break;
            }
            KillTime = KillTimes.at(KillCount);

            std::cout << "  Starting Checkpoint/Restart (CR) SixTrack run" << std::endl;

#ifndef WIN32
            pid_t SixTrackpid = fork();
            if(SixTrackpid == -1) {
                std::cerr << "ERROR: Could not fork to start SixTrack" << std::endl;
                return EXIT_FAILURE;
            }

            //Check fork() status
            if(SixTrackpid == 0) {
                //child, run sixtrack
                int execStatus = execl(argv[1], argv[1], (char*) 0);
                if(execStatus == -1) {
                    perror("ERROR: Could not execute SixTrack CR");
                }
            }
            else {
                std::cout << "  Will try to kill SixTrack after " << KillTime << " seconds" << std::endl;
                //Parent process
                pthread_t th_wait;  //This will call waitpid()
                pthread_t th_kill;  //This will do the kill()

                KillInfo th_wait_struct;
                th_wait_struct.SixPID = SixTrackpid;
                th_wait_struct.kTime = 0;
                th_wait_struct.RunStatus = false;    // True if SixTrack was finished
                th_wait_struct.CRKILLSWITCH = false; // True if SixTrack killed itself

                KillInfo th_kill_struct;
                th_kill_struct.SixPID = SixTrackpid;
                th_kill_struct.kTime = KillTime;
                th_kill_struct.RunStatus = false;
                th_kill_struct.CRKILLSWITCH = false;

                //The wait thread needs the pid + to return the status of the CR run
                // (killed or finished).
                int th_wait_ret = pthread_create(
                    &th_wait,
                    NULL,
                    pthread_wait_sixtrack,
                    (void*) &th_wait_struct
                );

                //The kill thread needs the pid to kill + the kill time for this cycle
                int th_kill_ret = pthread_create(
                    &th_kill,
                    NULL,
                    pthread_kill_sixtrack,
                    (void*) &th_kill_struct
                );

                pthread_join(th_kill, NULL);
                pthread_join(th_wait, NULL);

                if (th_wait_struct.RunStatus == true) {
                    std::cout << "  CR run finished! Will terminate the loop." << std::endl;
                    break;
                }
                if (th_wait_struct.CRKILLSWITCH == false) {
                    KillCount++;
                }
            } // END IF

#else
            //a DWORD is a 32bit unsigned int

            /*
              typedef struct _STARTUPINFO {
                DWORD  cb;
                LPTSTR lpReserved;
                LPTSTR lpDesktop;
                LPTSTR lpTitle;
                DWORD  dwX;
                DWORD  dwY;
                DWORD  dwXSize;
                DWORD  dwYSize;
                DWORD  dwXCountChars;
                DWORD  dwYCountChars;
                DWORD  dwFillAttribute;
                DWORD  dwFlags;
                WORD   wShowWindow;
                WORD   cbReserved2;
                LPBYTE lpReserved2;
                HANDLE hStdInput;
                HANDLE hStdOutput;
                HANDLE hStdError;
              } STARTUPINFO, *LPSTARTUPINFO;
            */
            STARTUPINFO SixTrack_si;

            /*
              typedef struct _PROCESS_INFORMATION {
                HANDLE hProcess;
                HANDLE hThread;
                DWORD  dwProcessId;
                DWORD  dwThreadId;
              } PROCESS_INFORMATION, *LPPROCESS_INFORMATION;
            */
            PROCESS_INFORMATION SixTrack_pi;

            //crashes happen if these are not first zero'ed
            ZeroMemory( &SixTrack_si, sizeof(SixTrack_si) );
            SixTrack_si.cb = sizeof(SixTrack_si);
            ZeroMemory( &SixTrack_pi, sizeof(SixTrack_pi) );

            /*
              BOOL WINAPI CreateProcess(
                _In_opt_    LPCTSTR               lpApplicationName,
                _Inout_opt_ LPTSTR                lpCommandLine,
                _In_opt_    LPSECURITY_ATTRIBUTES lpProcessAttributes,
                _In_opt_    LPSECURITY_ATTRIBUTES lpThreadAttributes,
                _In_        BOOL                  bInheritHandles,
                _In_        DWORD                 dwCreationFlags,
                _In_opt_    LPVOID                lpEnvironment,
                _In_opt_    LPCTSTR               lpCurrentDirectory,
                _In_        LPSTARTUPINFO         lpStartupInfo,
                _Out_       LPPROCESS_INFORMATION lpProcessInformation
                );
            */
            std::cout << "  CreateProcess()" << std::endl;
            //Start SixTrack running
            BOOL cp = CreateProcess(NULL, argv[1], NULL, NULL, FALSE, CREATE_NO_WINDOW, NULL, NULL, &SixTrack_si, &SixTrack_pi);

            if(!cp) {
                std::cerr << "ERROR: When starting SixTrack via CreateProcess()"  << std::endl;
                std::cerr << "Tried to run: " << argv[1] << " and error " << GetLastError() << std::endl;
            }
            else {
                //Our process should now be running
                std::cout << "  Will try to kill SixTrack after " << KillTime << " seconds" << std::endl;

                //Make 2 threads, one to just wait, one to kill
                KillInfo th_wait_struct;
                th_wait_struct.SixHANDLE    = SixTrack_pi.hProcess;
                th_wait_struct.kTime        = 0;
                th_wait_struct.RunStatus    = false; // True if SixTrack was finished
                th_wait_struct.CRKILLSWITCH = false; // True if SixTrack killed itself

                KillInfo th_kill_struct;
                th_kill_struct.SixHANDLE    = SixTrack_pi.hProcess;
                th_kill_struct.kTime        = KillTime;
                th_kill_struct.RunStatus    = false;
                th_kill_struct.CRKILLSWITCH = false;

                /*
                  HANDLE WINAPI CreateThread(
                    _In_opt_  LPSECURITY_ATTRIBUTES  lpThreadAttributes,
                    _In_      SIZE_T                 dwStackSize,
                    _In_      LPTHREAD_START_ROUTINE lpStartAddress,
                    _In_opt_  LPVOID                 lpParameter,
                    _In_      DWORD                  dwCreationFlags,
                    _Out_opt_ LPDWORD                lpThreadId
                    );
                */
                //Now make 2 threads
                HANDLE t_wait =  CreateThread(NULL, 0, winthread_wait_sixtrack, (LPVOID) &th_wait_struct, 0, NULL);
                HANDLE t_kill =  CreateThread(NULL, 0, winthread_kill_sixtrack, (LPVOID) &th_kill_struct, 0, NULL);

                WaitForSingleObject(t_wait, INFINITE);
                WaitForSingleObject(t_kill, INFINITE);
                CloseHandle(t_wait);
                CloseHandle(t_kill);

                CloseHandle(SixTrack_pi.hProcess);
                CloseHandle(SixTrack_pi.hThread);

                if(th_wait_struct.RunStatus == true) {
                    std::cout << "  CR run finished! Will terminate the loop." << std::endl;
                    break;
                }
                if (th_wait_struct.CRKILLSWITCH == false) {
                    KillCount++;
                }
            } // ENDIF
#endif

            //Wait a moment until the next run attempt is started?
            sleep(1);

        } // End WHILE(TRUE)

        std::cout << "  End CR run loop" << std::endl;
        std::cout << std::endl;

        // If a minimum number of restarts have been set, check this.
        if (atoi(argv[12]) > 0) {
            PrettyDivider("CHECKING Restarts");

            // Open metadata file
            std::ifstream MetaFile("sim_meta.dat", std::ifstream::in);
            std::string metaStr;
            bool foundCR_RestartCount = false;
            while(std::getline(MetaFile, metaStr)) {
                if (metaStr.substr(0,32) == "CR_RestartCount                 ") {
                    if(foundCR_RestartCount) {
                        std::cout << "  ERROR: Found two lines with 'CR_RestartCount' in 'sim_meta.dat'" << std::endl;
                        crRestartFail = true;
                    }
                    foundCR_RestartCount = true;

                    std::string metaSubStr = metaStr.substr(34,std::string::npos);
                    metaSubStr.erase( std::remove( metaSubStr.begin(), metaSubStr.end(), ' ' ),
                                      metaSubStr.end() );

                    if (metaSubStr.length() == 0) {
                        std::cout << "  ERROR: Reading 'sim_meta.dat'; "
                                  << "could not find the number of restarts in the line '"
                                  << metaStr << "'" << std::endl;
                        crRestartFail = true;
                    }

                    int numRestarts = atoi(metaSubStr.c_str());
                    std::stringstream numCRStatus;
                    if (numRestarts < atoi(argv[12])) {
                        numCRStatus << "CR_RestartCount: " << numRestarts << " < " << atoi(argv[12]);
                        crRestartFail = true;
                    } else {
                        numCRStatus << "CR_RestartCount: " << numRestarts << " >= " << atoi(argv[12]);
                    }
                    CheckPrint(numCRStatus.str(),!crRestartFail);
                }
            }
            if (not foundCR_RestartCount) {
                std::cout << "  Did not find 'CR_RestartCount' in 'sim_meta.dat'" << std::endl;
                crRestartFail = true;
            }
            MetaFile.close();

            std::cout << std::endl;

        }
    }

    //Normal run
    else {
        PrettyDivider("Normal Run");
#ifndef WIN32
        pid_t SixTrackpid = fork();
        if(SixTrackpid == -1) {
            std::cerr << "ERROR: Could not fork to start SixTrack" << std::endl;
            return EXIT_FAILURE;
        }

        //Check fork() status
        if(SixTrackpid == 0) {
            //In child process, run sixtrack

            // Redirect STDOUT to fort.6
            // Solution from http://stackoverflow.com/questions/20488574/output-redirection-using-fork-and-execl
            int fd_6 = open("fort.6",O_WRONLY|O_CREAT|O_TRUNC, S_IRUSR | S_IWUSR);
            if(fd_6 < 0) {
                perror("ERROR: Could not open file 'fort.6' for output");
                exit(EXIT_FAILURE);
            }
            dup2(fd_6,1);

            //Execute!
            int execStatus = execl(argv[1], argv[1], (char*) 0);
            if(execStatus == -1) {
                perror("ERROR: Could not execute SixTrack");
            }
        }
        else {
            //In main thread, wait()
            std::cout << "  Waiting for SixTrack to finish running. PID: " << SixTrackpid << std::endl;
            int waitpidStatus;
            waitpid(SixTrackpid, &waitpidStatus, WUNTRACED);
            std::cout << "  SixTrack finished running. Exit status: " << waitpidStatus << std::endl;
            std::cout << std::endl;

            if(waitpidStatus != 0) {
                //Print the last NLINES lines
                size_t NLINES = 10;
                PrettyDivider("Last Lines of Output");
                std::cout << std::endl;
                // Adapted from http://stackoverflow.com/questions/11876290/c-fastest-way-to-read-only-last-line-of-text-file
                std::ifstream fort6("fort.6");
                if (not fort6.is_open()) {
                    perror("ERROR: Could not open file 'fort.6' after SixTrack has finished");
                }
                else {
                    // go to one spot before the EOF
                    fort6.seekg(-1,fort6.end);
                    char ch;
                    size_t haslines = 0;
                    bool aborted_reading = false;
                    while (haslines<=NLINES) {
                        fort6.get(ch);
                        // In case the file contained only a single line
                        if((int)fort6.tellg() <= 1) {
                            fort6.seekg(0);
                            aborted_reading = true;
                            break;
                        }
                        // If the data was a newline
                        else if(ch == '\n') {
                            haslines += 1;
                            // Move to whatever is before the '\n'
                            fort6.seekg(-2,fort6.cur);
                        }
                        // If the data was neither a newline nor at the 0 byte
                        else {
                            // Move to the front of that data,
                            //  then to the front of the data before it
                            fort6.seekg(-2,fort6.cur);
                        }
                    }
                    if (not aborted_reading) {
                        // Skip to after the final random character and \n
                        fort6.seekg(+2,fort6.cur);
                    }
                    //OK, now we have moved NLINES up ahead; time to read and print again
                    std::string fort6_line;
                    while (std::getline(fort6,fort6_line)){
                        std::cout << fort6_line<<std::endl;
                    }

                    fort6.close();

                    PrettyDivider("End of Output");
                    std::cout << std::endl;
                }
            }
        }
#else
        int execStatus = _spawnl(_P_WAIT, argv[1], argv[1], (char*) 0);
        if(execStatus == -1) {
            perror("ERROR: Could not execute SixTrack");
        }
        std::cout << "  SixTrack finished running. Exit status: " << execStatus << std::endl;
        std::cout << std::endl;
#endif
    }

    /*
     * The next step is to check the output is valid.
     * A number of checks against reference output files can be made
     * Deal with fort.10 fort.90 and STF
     * Also look for extra checks
     */

    if(fort10) {
        PrettyDivider("CHECKING fort.10");
        fort10fail = CheckFort10(argv);
        std::cout << std::endl;
    }

    if(fort90) {
        PrettyDivider("CHECKING fort.90");
        fort90fail = CheckFort90(argv);
        std::cout << std::endl;
    }

    //Look at STF
    if(STF) {
        PrettyDivider("CHECKING singletrackfile.dat");
        STFfail = CheckSTF(argv);
        std::cout << std::endl;
    }

    //Look at extra_checks.txt
    ExtraChecksfail = PerformExtraChecks(extrachecks, argv[10],argv[11]);

    //Look at sixout.zip
#ifdef LIBARCHIVE
    if (sixoutzip) {
        PrettyDivider("CHECKING sixout.zip");

        //(Re-)create tmpdir folder
        struct stat st;
        int status;
        if (stat(tmpdir, &st) == 0) {
            if ( ! S_ISDIR(st.st_mode) ) {
                std::cerr << "ERROR: " << tmpdir << " exists, but is not a directory. Strange?!?" << std::endl;
                exit(EXIT_FAILURE);
            }
            std::cout << "  Folder '" << tmpdir << "' exits, deleting contents." << std::endl;

            // From http://stackoverflow.com/questions/11007494/how-to-delete-all-files-in-a-folder-but-not-delete-the-folder-using-nix-standar

            // These are data types defined in the "dirent" header
            DIR *theFolder = opendir(tmpdir);
            struct dirent *next_file;
            char filepath[256];

            while ( (next_file = readdir(theFolder)) != NULL ) {
                // build the path for each file in the folder
                int snprintf_err = snprintf(filepath, 256, "%s/%s", tmpdir, next_file->d_name);
                if (snprintf_err >= 256 || snprintf_err < 0) {
                    std::cout << "ERROR: in snprintf while building the path for deleting a file" << std::endl;
                    exit(EXIT_FAILURE);
                }
                remove(filepath);
            }
            closedir(theFolder);
            std::cout << std::endl;
        }
        else {
            std::cout << "  Creating folder '" << tmpdir << "' ..." << std::endl;
#if defined(_WIN32)
            status = not CreateDirectory(tmpdir,NULL);
#else
            status = mkdir(tmpdir,S_IRWXU);
#endif
            if (status) {
                std::cerr << "ERROR: Something went wrong when creating '" << tmpdir << "'. Sorry!" << std::endl;
                exit(EXIT_FAILURE);
            }
            std::cout << std::endl;
        }

        //List content
        //list_archive(sixoutzip_fname);
        const int archive_nfiles_max = 256;
        int archive_nfiles = archive_nfiles_max;
        const int archive_buffsize = 100;
        char** archive_buff = new char*[archive_nfiles];
        for (int i = 0; i< archive_nfiles_max; i++) {
            archive_buff[i] = new char[archive_buffsize];
        }
        std::cout << "  Calling list_archive_get..." << std::endl;
        list_archive_get(sixoutzip_fname,archive_buff,&archive_nfiles,archive_buffsize);

        std::cout << "  Got " << archive_nfiles << " files:" << std::endl;
        for (int i = 0; i< archive_nfiles; i++) {
            std::cout << "  File #" << i << ": '" << archive_buff[i] << "'" << std::endl;
        }
        if (archive_nfiles != sixoutzip) {
            std::cout << "  WARNING: Expected " << sixoutzip << " files in '"
                      << sixoutzip_fname << "', got " << archive_nfiles << std::endl;
            sixoutzipfail = true;
        }
        else {
            std::cout << "  The number of files in the archive '"
                      << sixoutzip_fname << "' MATCHES the expected number "
                      << archive_nfiles << std::endl;
        }

        //Unzip!
        std::cout << "  Calling read_archive..." << std::endl;
        read_archive(sixoutzip_fname,tmpdir);

        for (int i = 0; i< archive_nfiles; i++) {
            char* FileNameZip = new char[strlen(tmpdir)+1+archive_buffsize];

            // Insert the right path separator
            int snprintf_err = 0;
#ifdef WIN32
            snprintf_err = snprintf(FileNameZip,archive_buffsize+1+strlen(tmpdir),"%s\\%s",tmpdir,archive_buff[i]);
#else
            snprintf_err = snprintf(FileNameZip,archive_buffsize+1+strlen(tmpdir),"%s/%s",tmpdir,archive_buff[i]);
#endif
            if (snprintf_err >= (archive_buffsize+1+strlen(tmpdir)) || snprintf_err < 0) {
                std::cerr << "ERROR: In snprintf, while building the path for unzipped file."
                          << std::endl;
                std::cerr << "snprintf_err = " << snprintf_err
                          << ", buffsize=" << (archive_buffsize+1+strlen(tmpdir)) << std::endl;
                exit(EXIT_FAILURE);
            }

#ifdef WIN32
            //Strip out \r characters from windows new lines
            size_t CRcount = StripCR(FileNameZip);
            std::cout << "  Removed " << CRcount << " windows \\r entries from '"
                      << FileNameZip << "'." << std::endl;
#endif

            bool ThisTest = !FileComparison(FileNameZip, std::string(archive_buff[i]) + ".canonical");
            if(ThisTest) {
                std::cerr << "  WARNING: Test of zipped file '" << FileNameZip << "' failed!" << std::endl;
                sixoutzipfail = true;
            }
            else {
                std::cout << "  Test of zipped file '" << FileNameZip << "' MATCHES" << std::endl;
            }
        }

        //Cleanup memory
        for (int i = 0; i< archive_nfiles_max; i++) {
            delete[] archive_buff[i];
        }
        delete[] archive_buff;
        std::cout << std::endl;
    }
#endif
    PrettyDivider("CHECKS SUMMARY");
    if(fort10)      CheckPrint("fort.10",!fort10fail);
    if(fort90)      CheckPrint("fort.90",!fort90fail);
    if(STF)         CheckPrint("singletrackfile.dat",!STFfail);
    if(extrachecks) CheckPrint("Extra Checks",!ExtraChecksfail);
#ifdef LIBARCHIVE
    if(sixoutzip)   CheckPrint(sixoutzip_fname,!sixoutzipfail);
#endif
    if(CRon && atoi(argv[12]) > 0) {
        CheckPrint("CR Number of Restarts",!crRestartFail);
    }
    std::cout << std::endl;

    PrettyDivider("EXIT");
    std::cout << std::endl;

    //or together any fail bits.
    //If all tests pass this will return 0 (good)
    //if not we get something else out (bad)
    return (fort10fail || fort90fail || STFfail || ExtraChecksfail || sixoutzipfail || crRestartFail);
}

/**
* Copies one input file to an output file.
* @parm InputFileName The name of the input file
* @parm OutputFileName The name of the output file
* @return A bool set to true if the operation succeeded
*/
bool CopyFile(std::string InputFileName, std::string OutputFileName) {
    std::ifstream Input(InputFileName.c_str(), std::ios::binary);
    std::ofstream Output(OutputFileName.c_str(), std::ios::binary);

    if(!Input.good()) {
        std::cerr << "  ERROR: Could not open " << InputFileName << std::endl;
        return false;
    }

    if(!Input.good()) {
        std::cerr << "  ERROR: Could not open " << OutputFileName << std::endl;
        return false;
    }

    Output << Input.rdbuf();
    return true;
}

void CheckPrint(std::string theText, bool checkResult) {
    int nMax = 80;
    int nDot = nMax-theText.length()-14;
    if(checkResult) {
        std::cout << "  " << theText << " ";
        for(int i=0; i<nDot; i++) {
            std::cout << ".";
        }
        std::cout << "   Passed" << std::endl;
    } else {
        std::cerr << "  " << theText << " ";
        for(int i=0; i<nDot; i++) {
            std::cerr << ".";
        }
        std::cerr << "***Failed" << std::endl;
    }
}

void PrettyDivider(std::string theText) {
    int nMax = 80;
    int nSt  = (nMax-2-theText.length())/2;
    int nEnd = nMax-nSt-2-theText.length();
    if(theText.length() == 0) {
        for(int i=0; i<nMax; i++) {
            std::cout << "*";
        }
    } else {
        for(int i=0; i<nSt; i++) {
            std::cout << "*";
        }
        std::cout << " " << theText << " ";
        for(int i=0; i<nEnd; i++) {
            std::cout << "*";
        }
    }
    std::cout << std::endl;
}

/**
* Checks if two fort.10 files are equivalent
* @return true if there is a problem, false if there is not
*/
bool CheckFort10(char* argv[]) {
    //return false if all is good, true if anything else happens

    //In this folder we should have fort.10 and a fort.10.canonical
    bool input1 = CopyFile("fort.10", "fort.20");
    bool input2 = CopyFile("fort.10.canonical", "fort.21");

    if(!input1 || !input2) {
        std::cerr << "  WARNING: Could not perform fort.10 comparison" << std::endl;
        return true;
    }

#ifndef WIN32
    //Now again we fork() and exec()
    pid_t CheckFort10pid = fork();
    if(CheckFort10pid == -1) {
        std::cerr << "  ERROR: Could not fork to start checkf10" << std::endl;
        return true;
    }

    //Check fork() status
    if(CheckFort10pid == 0) {
        //child, run checkf10
        int status = execl(argv[2], "checkf10", (char*) 0);
        if(status == -1) {
            perror("ERROR: Could not execute checkf10");
            exit(EXIT_FAILURE);
        }
    }
    else {
        //main thread, wait()
        int waitpidStatus;
        waitpid(CheckFort10pid, &waitpidStatus, WUNTRACED);
        std::cout << "  checkf10 finished running and returned: " << waitpidStatus << std::endl;
        return waitpidStatus;
    }

    //If we hit here, there is a problem!
    return true;
#else
    int status = _spawnl(_P_WAIT, argv[2], "checkf10", (char*) 0);
    if(status == -1){
        perror("ERROR: Could not execute checkf10");
    }
    std::cout << "  checkf10 finished running and returned: " << status << std::endl;

    return status;
#endif
}

/**
* Checks if two fort.90 files are equivalent
* @return true if there is a problem, false if there is not
*/
bool CheckFort90(char* argv[]) {
    //return false if all is good, true if anything else happens
    //First we call read90 on each file.
    //Then we must do a binary comparison

#ifndef WIN32
    //Now again we fork() and exec()
    //Do this for the first file
    pid_t CheckFort90pid = fork();
    if(CheckFort90pid == -1) {
        std::cerr << "  ERROR: Could not fork to start read90" << std::endl;
        return false;
    }

    //Check fork() status
    if(CheckFort90pid == 0){
        //child, run read90
        int status = execl(argv[3], "read90",
                           "--fname", "fort.90",
                           "--ofname", "fort.90.out", (char*) 0);
        if(status == -1) {
            perror("ERROR: Could not execute read90");
            exit(EXIT_FAILURE);
        }
    }
    else {
        //main thread, wait()
        int waitpidStatus;
        waitpid(CheckFort90pid, &waitpidStatus, WUNTRACED);
        std::cout << "  read90 finished running on fort.90: " << waitpidStatus << std::endl;
        if(waitpidStatus != 0) {
            std::cerr << "  ERROR: Problem running read90" << std::endl;
            return waitpidStatus;
        }
    }

    //Do this for the second file
    CheckFort90pid = fork();
    if(CheckFort90pid == -1) {
        std::cerr << "  ERROR: Could not fork to start read90" << std::endl;
        return false;
    }

    //Check fork() status
    if(CheckFort90pid == 0) {
        //child, run read90
        int status = execl(argv[3], "read90",
                           "--fname", "fort.90.canonical",
                           "--ofname", "fort.90.canonical.out", (char*) 0);
        if(status == -1) {
            perror("ERROR: Could not execute read90");
            exit(EXIT_FAILURE);
        }
    }
    else {
        //main thread, wait()
        int waitpidStatus;
        waitpid(CheckFort90pid, &waitpidStatus, WUNTRACED);
        std::cout << "  read90 finished running on fort.90.canonical: " << waitpidStatus << std::endl;
        if(waitpidStatus != 0) {
            std::cerr << "  ERROR: Problem running read90" << std::endl;
            return waitpidStatus;
        }
    }
#else
    //Windows
    //First file
    int status1 = _spawnl(_P_WAIT, argv[3], "read90",
                          "--fname", "fort.90",
                          "--ofname", "fort.90.out", (char*) 0);
    if(status1 == -1){
        perror("ERROR: Could not execute read90");
    }
    std::cout << "  read90 finished running on fort.90: " << status1 << std::endl;

    //Second file (canonical)
    int status2 = _spawnl(_P_WAIT, argv[3], "read90",
                          "--fname", "fort.90.canonical",
                          "--ofname", "fort.90.canonical.out", (char*) 0);
    if(status2 == -1) {
        perror("ERROR: Could not execute read90");
    }
    std::cout << "  read90 finished running on fort.90.canonical: " << status2 << std::endl;
#endif
    return !FileComparison("fort.90.out", "fort.90.canonical.out");

}

/**
* Checks if two singletrackfile.dat files are equivalent
* @return true if there is a problem, false if there is not
*/
bool CheckSTF(char* argv[]) {
    //return false if all is good, true if anything else happens
    //First we call read90 on each file.
    //Then we must do a binary comparison
#ifndef WIN32
    //Now again we fork() and exec()
    //Do this for the first file
    pid_t CheckSTFpid = fork();
    if(CheckSTFpid == -1) {
        std::cerr << "  ERROR: Could not fork to start read90" << std::endl;
        return false;
    }

    //Check fork() status
    if(CheckSTFpid == 0) {
        //child, run read90
        int status = execl(argv[3], "read90", "--STF",
                           "--fname", "singletrackfile.dat",
                           "--ofname", "singletrackfile.dat.out", (char*) 0);
        if(status == -1) {
            perror("ERROR: Could not execute read90");
            exit(EXIT_FAILURE);
        }
    }
    else {
        //main thread, wait()
        int waitpidStatus;
        waitpid(CheckSTFpid, &waitpidStatus, WUNTRACED);
        std::cout << "  read90 finished running on singletrackfile.dat: "
                  << waitpidStatus << std::endl;
        if(waitpidStatus != 0) {
            std::cerr << "  ERROR: Problem running read90" << std::endl;
            return waitpidStatus;
        }
    }

    //Do this for the second file
    CheckSTFpid = fork();
    if(CheckSTFpid == -1) {
        std::cerr << "  ERROR: Could not fork to start read90" << std::endl;
        return false;
    }

    //Check fork() status
    if(CheckSTFpid == 0) {
        //child, run read90
        int status = execl(argv[3], "read90", "--STF",
                           "--fname", "singletrackfile.dat.canonical",
                           "--ofname", "singletrackfile.dat.canonical.out", (char*) 0);
        if(status == -1) {
            perror("ERROR: Could not execute read90");
            exit(EXIT_FAILURE);
        }
    }
    else {
        //main thread, wait()
        int waitpidStatus;
        waitpid(CheckSTFpid, &waitpidStatus, WUNTRACED);
        std::cout << "  read90 finished running on singletrackfile.dat.canonical: "
                  << waitpidStatus << std::endl;
        if(waitpidStatus != 0) {
            std::cerr << "  ERROR: Problem running read90" << std::endl;
            return waitpidStatus;
        }
    }
#else
    //Windows
    //First file
    int status1 = _spawnl(_P_WAIT, argv[3], "read90", "--STF",
                          "--fname", "singletrackfile.dat",
                          "--ofname", "singletrackfile.dat.out", (char*) 0);
    if(status1 == -1) {
        perror("ERROR: Could not execute read90");
    }
    std::cout << "  read90 finished running on singletrackfile.dat: "
              << status1 << std::endl;

    //Second file (canonical)
    int status2 = _spawnl(_P_WAIT, argv[3], "read90", "--STF",
                          "--fname", "singletrackfile.dat.canonical",
                          "--ofname", "singletrackfile.dat.canonical.out", (char*) 0);
    if(status2 == -1) {
        perror("ERROR: Could not execute read90");
    }
    std::cout << "  read90 finished running on singletrackfile.dat.canonical: "
              << status2 << std::endl;
#endif
    return !FileComparison("singletrackfile.dat.out", "singletrackfile.dat.canonical.out");

}
/**
* File comparison
*/
bool FileComparison(std::string FileName1, std::string FileName2) {
    std::ifstream f1(FileName1.c_str(), std::ios::binary);
    std::ifstream f2(FileName2.c_str(), std::ios::binary);

    size_t length1 = 0;
    size_t length2 = 0;
    if(f1.good()) {
        f1.seekg(0,f1.end);
        length1 = f1.tellg();
    }
    else {
        std::cout << "  Could not open " << FileName2 << std::endl;
        return false;
    }

    if(f2.good()) {
        f2.seekg(0,f2.end);
        length2 = f2.tellg();
    }
    else {
        std::cout << "  Could not open " << FileName2 << std::endl;
        return false;
    }

    std::cout << "  Size " << FileName1 << ": " << length1 << std::endl;
    std::cout << "  Size " << FileName2 << ": " << length2 << std::endl;

    if(length1 != length2) {
        return false;
    }

    //Same size files, now to check the contents
    const size_t BufferSize = 1024;
    char f1Buffer[BufferSize];
    char f2Buffer[BufferSize];

    f1.seekg(0);
    f2.seekg(0);

    size_t position = 0;

    while(position <= length1) {
        size_t rsize = std::min(length1, BufferSize);
        f1.read(f1Buffer,rsize);
        f2.read(f2Buffer,rsize);
        //std::cerr << FileName1 << " is different from " << FileName2
        //          << " at " << length2-length1 << " - comparison: " << comparison << std::endl;

        int comparison = memcmp(f1Buffer, f2Buffer, rsize);
        if(comparison != 0) {
            std::cerr << "  " << FileName1 << " is different from " << FileName2
                      << " at " << position << " - memcmp() comparison: " << comparison << std::endl;
            return false;
        }

        position+= rsize;
    }
    return true;
}

bool PerformExtraChecks(bool &extrachecks, char* convert_dump_bin, char* dump_bin_files) {

    PrettyDivider("EXTRA CHECKS");

    //Split the dump_bin_files
    std::list<std::string> dump_bin_files_list;
    std::string dump_bin_files_input = dump_bin_files;
    size_t dump_bin_files_pos = 0;
    size_t dump_bin_files_oldpos = 0;
    while(dump_bin_files_pos != std::string::npos) {
        dump_bin_files_pos = dump_bin_files_input.find(",", dump_bin_files_oldpos);
        dump_bin_files_list.push_back( dump_bin_files_input.substr(dump_bin_files_oldpos,
                                                                   dump_bin_files_pos
                                                                    - dump_bin_files_oldpos)
                                       );
        dump_bin_files_oldpos=dump_bin_files_pos+1;
        std::cout << "  Found binary dump file = '" << dump_bin_files_list.back() << "'" << std::endl;
    }
    if (dump_bin_files_list.size()==1 and dump_bin_files_list.front()=="NONE") {
        dump_bin_files_list.clear();
    }

    bool AllTests = false;
    std::ifstream extra_checks_in("extra_checks.txt");
    if(extra_checks_in.good()) {
        extrachecks=true;
        std::cout << "  Opened extra_checks.txt" << std::endl;
        //Format should be some file to check followed by a command
        while(extra_checks_in.good()) {
            std::string StringBuffer;
            std::string FileName;

            extra_checks_in >> FileName;
            getline(extra_checks_in, StringBuffer);
            if(FileName != "") {
                std::cout << "  Performing extra checks on '" << FileName << "'" << std::endl;
                bool convertThis = false;
                auto dump_bin_files_iterator = std::find(dump_bin_files_list.begin(),dump_bin_files_list.end(),FileName);
                if(dump_bin_files_iterator != dump_bin_files_list.end()) {
                    convertThis = true;
                    std::cout << "  This is a binary format 3 DUMP, must convert!" << std::endl;
                    std::cout << "  Calling '"<<convert_dump_bin<<"'"<<std::endl;
#ifndef WIN32
                    //Now again we fork() and exec()
                    //Do this for the first file
                    pid_t ReadDump3pid = fork();
                    if(ReadDump3pid == -1) {
                        perror("ERROR: Could not fork to start readDump3");
                        exit(EXIT_FAILURE);
                    }

                    //Check fork() status
                    if(ReadDump3pid == 0) {
                        //child, run readDump3
                        int status = execl(convert_dump_bin, "readDump3",
                                           FileName.c_str(),
                                           (FileName+std::string(".converted")).c_str(),
                                           (char*) 0);
                        if(status == -1) {
                            perror("ERROR: Could not execute readDump3");
                            exit(EXIT_FAILURE);
                        }
                    }
                    else {
                        //main thread, wait()
                        int waitpidStatus;
                        waitpid(ReadDump3pid, &waitpidStatus, WUNTRACED);
                        std::cout << "  readDump3 finished running on '"<< FileName << "': "
                                  << waitpidStatus << std::endl;
                        if(waitpidStatus != 0) {
                            std::cerr << "ERROR: Problem running readDump3" << std::endl;
                            exit(EXIT_FAILURE);
                        }
                    }

                    //Do this for the second file
                    ReadDump3pid = fork();
                    if(ReadDump3pid == -1) {
                        perror("ERROR: Could not fork to start readDump3");
                        exit(EXIT_FAILURE);
                    }

                    //Check fork() status
                    if(ReadDump3pid == 0) {
                        //child, run readDump3
                        int status = execl(convert_dump_bin, "readDump3",
                                           (FileName.c_str()+std::string(".canonical")).c_str(),
                                           (FileName+std::string(".converted")
                                            + std::string(".canonical")).c_str(),
                                           (char*) 0);
                        if(status == -1) {
                            perror("ERROR: Could not execute readDump3");
                            exit(EXIT_FAILURE);
                        }
                    }
                    else {
                        //main thread, wait()
                        int waitpidStatus;
                        waitpid(ReadDump3pid, &waitpidStatus, WUNTRACED);
                        std::cout << "  readDump3 finished running on '"
                                  << FileName.c_str()+std::string(".canonical") << "': "
                                  << waitpidStatus << std::endl;
                        if(waitpidStatus != 0) {
                            std::cerr << "ERROR: Problem running readDump3" << std::endl;
                            exit(EXIT_FAILURE);
                        }
                    }
#else
                    //Windows
                    //First file
                    int status1 = _spawnl(_P_WAIT, convert_dump_bin, "readDump3",
                                          FileName.c_str(),
                                          (FileName+std::string(".converted")).c_str(),
                                          (char*) 0);
                    if(status1 == -1) {
                        perror("ERROR: Could not execute readDump3");
                        exit(EXIT_FAILURE);
                    }
                    std::cout << "  readDump3 finished running on '"
                              << FileName << "': "
                              << status1 << std::endl;

                    //Second file (canonical)
                    int status2 = _spawnl(_P_WAIT, convert_dump_bin, "read90",
                                          (FileName.c_str()+std::string(".canonical")).c_str(),
                                          (FileName+std::string(".converted")
                                            + std::string(".canonical")).c_str(),
                                          (char*) 0
                                          );
                    if(status2 == -1) {
                        perror("ERROR: Could not execute readDump3");
                    }
                    std::cout << "  readDump3 finished running on '"
                              << FileName.c_str() + std::string(".canonical")
                              << "': " << status2 << std::endl;
#endif
                    // Update the filename
                    FileName = FileName+std::string(".converted");
                } // Done converting to ASCII...
#ifdef WIN32
                //Strip out \r characters from windows new lines
                size_t CRcount = StripCR(FileName);
                std::cout << "  Removed " << CRcount << " windows \\r entries." << std::endl;
                if (convertThis) {
                    size_t CRcount = StripCR(FileName + std::string(".canonical"));
                    std::cout << "  Removed " << CRcount << " windows \\r entries." << std::endl;
                }
#endif
                bool ThisTest = !FileComparison(FileName, FileName + ".canonical");
                CheckPrint(FileName,!ThisTest);
                if(ThisTest) AllTests = true;
            }
        }
    }
    else {
        std::cout << "  Could not open extra_checks.txt" << std::endl;
        extrachecks=false;
    }

    std::cout << std::endl;
    return AllTests;
}

std::vector<int> ParseKillTimes(char* in) {
    std::vector<int> KillTimes;
    std::string input = in;
    size_t pos = 0;
    size_t oldpos = 0;
    while(pos != std::string::npos) {
        pos = input.find(",", oldpos);
        int number = atoi(input.substr(oldpos,pos-oldpos).c_str());
        oldpos=pos+1;

        KillTimes.push_back(number);
    }

    PrettyDivider("Checkpoint/Restart");
    std::cout << "  Will try and kill CR run " << KillTimes.size() << " times." << std::endl;
    std::cout << "  Killing after ";
    for(int k = 0; k < KillTimes.size(); k++) {
        std::cout << KillTimes.at(k) << ", ";
    }
    std::cout << " seconds." << std::endl;
    return KillTimes;
}

#ifndef WIN32
void *pthread_wait_sixtrack(void* InputStruct) {
    //Grab the structure
    KillInfo* ThreadStruct = (KillInfo*)InputStruct;

    //Extract the PID from the structure
    pid_t sixpid = ThreadStruct->SixPID;

    int waitpidStatus;
    //Wait for the thread to either end or to be killed.
    waitpid(sixpid, &waitpidStatus, WUNTRACED);

    //If the thread exited (and was not killed), then we can exit this run.
    if(waitpidStatus == 0) {
        //Check if it was stopped by CRKILLSWITCH; if so remove the file
        int unlink_status = unlink("crrestartme.tmp");
        if (unlink_status != 0) {
            std::cout << "  SixTrack CR exited okay: " << waitpidStatus << std::endl;
            ThreadStruct->RunStatus    = true;
            ThreadStruct->CRKILLSWITCH = false;
        }
        else {
            std::cout << "  SixTrack CR was stopped by CRKILLSWITCH; 'crrestartme.tmp' was deleted."
                      << std::endl;
            ThreadStruct->RunStatus    = false;
            ThreadStruct->CRKILLSWITCH = true;
        }
    }
    else {
        std::cout << "  SixTrack CR was killed: " << waitpidStatus << std::endl;
        ThreadStruct->RunStatus    = false;
        ThreadStruct->CRKILLSWITCH = false;
    }

    pthread_exit(NULL);
}

void *pthread_kill_sixtrack(void* InputStruct) {
    //Grab the structure
    KillInfo* ThreadStruct = (KillInfo*)InputStruct;

    //Extract the kill time from the structure
    int KillTime = ThreadStruct->kTime;

    //Extract the PID from the structure
    pid_t sixpid = ThreadStruct->SixPID;
    bool ArmKill=true;

    //Countdown clock
    for(int tt=0; tt < KillTime; tt++) {
        sleep(1);
        //std::cout << "At " << tt+1 << " of " << KillTime << " Testing pid " << sixpid << ": ";

        int res = kill(sixpid, 0); //Signal 0 means "check if we could kill this if we wanted to".
        //std::cout << "Child exec() kill 0 check result: " <<  res << std::endl;
        if(res != 0) {

            std::cout << "  Cannot kill SixTrack process, most likely it has finished." << std::endl;

            std::cout << "  Error message from kill(): '" << strerror(errno) << "'" << std::endl;

            //No longer running, jump out;
            ArmKill=false;
            tt=KillTime;
        }
    }

    // Do the actual kill
    if(ArmKill == true) {
        //Try and kill
        std::cout << "  Kill thread - killing {IF}: " << sixpid << std::endl;
        int res = kill(sixpid, SIGKILL);
        std::cout << "  Kill thread - kill() result: " << res << std::endl;
    }

    pthread_exit(NULL);
}
#else
DWORD winthread_wait_sixtrack(LPVOID InputStruct) {
    //Grab the structure
    KillInfo* ThreadStruct = (KillInfo*)InputStruct;

    //Extract the HANDLE from the structure
    HANDLE sixHANDLE = ThreadStruct->SixHANDLE;

    //Wait for the thread to either end or to be killed.
    WaitForSingleObject(sixHANDLE, INFINITE);

    //If the thread exited (and was not killed), then we can exit this run.
    /*
      BOOL WINAPI GetExitCodeProcess(
        _In_  HANDLE  hProcess,
        _Out_ LPDWORD lpExitCode
      );
    */
    DWORD excode = 2;
    LPDWORD excode_ptr = &excode;
    BOOL ecode = GetExitCodeProcess(sixHANDLE, excode_ptr);
    if(!ecode) {
        std::cerr << "  GetExitCodeProcess() failed on SixTrack CR run!: " << GetLastError() << std::endl;
    }
    else {
        if(excode == 0) {
            //Check if it was stopped by CRKILLSWITCH; if so remove the file
            int unlink_status = unlink("crrestartme.tmp");
            if (unlink_status != 0) {
                std::cout << "  SixTrack CR exited okay: " << excode << std::endl;
                ThreadStruct->RunStatus    = true;
                ThreadStruct->CRKILLSWITCH = false;
            }
            else {
                std::cout << "  SixTrack CR was stopped by CRKILLSWITCH; 'crrestartme.tmp' was deleted."
                          << std::endl;
                ThreadStruct->RunStatus    = false;
                ThreadStruct->CRKILLSWITCH = true;
            }
        }
        else {
            std::cout << "  SixTrack CR was killed: " << excode << std::endl;
            ThreadStruct->RunStatus    = false;
            ThreadStruct->CRKILLSWITCH = false;
        }
    }

    return 0;
}

DWORD winthread_kill_sixtrack(LPVOID InputStruct) {
    //Grab the structure
    KillInfo* ThreadStruct = (KillInfo*)InputStruct;

    //Extract the kill time from the structure
    int KillTime = ThreadStruct->kTime;

    //Extract the HANDLE from the structure
    HANDLE sixHANDLE = ThreadStruct->SixHANDLE;
    bool ArmKill=true;

    for(int tt=0; tt < KillTime; tt++) {
        Sleep(1000);
        //std::cout << "At " << tt+1 << " of " << KillTime << " Testing pid " << sixpid << ": ";
        DWORD excode = 2;
        LPDWORD excode_ptr = &excode;
        BOOL ecode = GetExitCodeProcess(sixHANDLE, excode_ptr);
        if(!ecode) {
            std::cerr << "  GetExitCodeProcess() failed on SixTrack CR run!: "
                      << GetLastError() << std::endl;
        }
        else {
            if(excode == 0) {
                std::cout << "  SixTrack CR exited okay: " << excode << std::endl;
                //No longer running, jump out;
                ArmKill=false;
                tt=KillTime;
            }
        }
    }
    if(ArmKill == true) {
        //Try and kill
        std::cout << "  Kill thread - calling TerminateProcess()" << std::endl;
        /*
          BOOL WINAPI TerminateProcess(
            _In_ HANDLE hProcess,
            _In_ UINT   uExitCode
            );
        */
        BOOL term = TerminateProcess(sixHANDLE, 9);
        if(!term) {
            std::cerr << "  Failed to TerminateProcess() on SixTrack CR run!: "
                      << GetLastError() << std::endl;
        }
        std::cout << "  Kill thread - TerminateProcess() result: " << term << std::endl;
    }

    return 0;
}
#endif

/**
* Deletes any checkpoint files that are appended to from previous CR runs.
*/
void UnlinkCRFiles() {

    std::vector<std::string> unlinkFiles;

    unlinkFiles.push_back("fort.6");
    unlinkFiles.push_back("fort.10");
    unlinkFiles.push_back("fort.90");
    unlinkFiles.push_back("cr_status.log");
    unlinkFiles.push_back("crpoint_pri.bin");
    unlinkFiles.push_back("crpoint_sec.bin");
    unlinkFiles.push_back("crrestartme.tmp");
    unlinkFiles.push_back("crkillswitch.tmp");

    for(auto fname : unlinkFiles) {
        std::cout << "  Deleting old '" << fname << "'";
        if(unlink(fname.c_str()) != 0) {
            std::cout << " - File not found";
        } else {
            std::cout << " - Done";
        }
        std::cout << std::endl;
    }
}

size_t StripCR(std::string FileName) {
    std::ifstream InputFile(FileName.c_str(), std::ios::binary);
    if(InputFile.good()) {
        size_t count = 0;
        char cbuffer;
        std::streambuf* sbuffer = InputFile.rdbuf();
        std::stringstream stream;

        while(sbuffer->sgetc() != EOF) {
            if( (cbuffer = sbuffer->sbumpc()) ) {
                if(cbuffer == '\r') {
                    //Increment the count and do nothing
                    count++;
                }
                else {
                    stream << cbuffer;
                }
            }
        }

        //Now close the input file
        InputFile.close();

        //write!
        std::ofstream OutputFile(FileName.c_str(), std::ios::binary);
        OutputFile << stream.str();
        OutputFile.close();

        return count;
    }

    return 0;
}
