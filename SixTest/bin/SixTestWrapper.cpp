#include <cstdlib>
#include <iostream>
#include <fstream>

#include <string>
#include <vector>

#include <unistd.h>
#include <signal.h>

//waitpid()
#include <sys/types.h>
#include <sys/wait.h>

//to get errors from POSIX calls
#include <errno.h>
#include <stdio.h>

//memcmp
#include <string.h>

void KillSixTrack(size_t KillTime, pid_t sixpid);

bool CopyFile(std::string InputFileName, std::string OutputFileName);
bool FileComparison(std::string f1, std::string f2);

bool CheckFort10();
bool CheckFort90();
bool CheckSTF();
bool PerformExtraChecks();
std::vector<int> ParseKillTimes(char*);

/**
* SixTrack testing wrapper
* This runs a given SixTrack binary (given as the first argument)
* It then performs checks on the output
* Currently we check fort.10 and fort.90
* Arguments:
* 1: Sixtrack binary to run
* 2: bool to check fort.10
* 3: bool to check fort.90
* 4: bool for STF enabled
* 5: CR enabled
* 6: CR kill time
*/
int main(int argc, char* argv[])
{
	//First check we have the correct number of arguments
	if(argc != 7)
	{
		std::cout << argv[0] << " called with the incorrect number of arguments, should be 6, but was called with " << argc - 1 << " arguments" << std::endl;
		return EXIT_FAILURE;
	}

	//Want our first argument to be the sixtrack binary name, assume this is run in the test folder.
	std::cout << "Called with " << argc << " arguments" << std::endl;
	for(int n=0; n < argc; n++)
	{
		std::cout << "Argument: " << n << " is " << argv[n] << std::endl;
	}

	//Set defaults:
	//Set to true if we have a CR build
	bool CR = false;

	//How long to wait in seconds before killing the CR run
	std::vector<int> KillTimes;
	int KillTime = 0;

	bool fort6 = false;
	bool fort10 = false;
	bool fort90 = false;
	bool STF = false;

	bool fort10fail = false;
	bool fort90fail = false;
	bool STFfail = false;
	bool ExtraChecksfail = false;

	if(atof(argv[2]) != 0)
	{
		fort10 = true;
	}

	if(atof(argv[3]) != 0)
	{
		fort90 = true;
	}

	if(atof(argv[4]) != 0)
	{
		STF = true;
	}

	if(atof(argv[5]) != 0)
	{
		CR = true;
		KillTimes = ParseKillTimes(argv[6]);
	}

	/**
	* First step is to handle all the running of sixtrack.
	* This include Checkpoint/resume builds where the run will be killed after a set period of time and then restarted.
	* In each case, to run sixtrack we call fork and let the child process do the exec()
	* The main thread then calls waitpid() for sixtrack to finish executing.
	* For CR we wait then send a kill signal after a specified number of seconds
	*/
	if(CR)
	{
		for(size_t KillCount=0; KillCount < KillTimes.size(); KillCount++)
		{
			KillTime = KillTimes.at(KillCount);

			std::cout << "Starting Checkpoint/Resume (CR) SixTrack run" << std::endl;
			pid_t SixTrackpid = fork();
			if(SixTrackpid == -1)
			{
				std::cerr << "ERROR: Could not fork to start SixTrack" << std::endl;
				return EXIT_FAILURE;
			}

			//Check fork() status
			if(SixTrackpid == 0)
			{
				//child, run sixtrack
				int execStatus = execl(argv[1], argv[1], (char*) 0);
				if(execStatus == -1)
				{
					perror("ERROR - could not execute SixTrack");
				}
			}
			else
			{
				std::cout << "Will kill SixTrack after " << KillTime << " seconds" << std::endl;
				KillSixTrack(KillTime, SixTrackpid);
			}
		}
	}

	//Normal run
	else
	{
		pid_t SixTrackpid = fork();
		if(SixTrackpid == -1)
		{
			std::cerr << "ERROR: Could not fork to start SixTrack" << std::endl;
			return EXIT_FAILURE;
		}

		//Check fork() status
		if(SixTrackpid == 0)
		{
			//child, run sixtrack
			int execStatus = execl(argv[1], argv[1], (char*) 0);
			if(execStatus == -1)
			{
				perror("ERROR - could not execute SixTrack");
			}
		}
		else
		{
			//main thread, wait()
			std::cout << "Waiting for SixTrack to finish running - pid: " << SixTrackpid << std::endl;
			int waitpidStatus;
			waitpid(SixTrackpid, &waitpidStatus, WUNTRACED);
			std::cout << "SixTrack finished running: " << waitpidStatus << std::endl;
		}
	}
	/*
	* The next step is to check the output is valid.
	* A number of checks against reference output files can be made
	* Deal with fort.10 fort.90 and STF
	* Also look for extra checks
	*/
	std::cout << "Will now check output files" << std::endl;
	if(fort10)
	{
		std::cout << "------------------------------ Checking fort.10 ------------------------------" << std::endl;
		fort10fail = CheckFort10();
		std::cout << "---------------------------- End checking fort.10 ----------------------------" << std::endl;
		if(fort10fail)
		{
			std::cerr << "WARNING: fort.10 files do NOT match" << std::endl;
		}
		else
		{
			std::cout << "fort.10: PASS" << std::endl;
		}
	}

	if(fort90)
	{
		std::cout << "------------------------------ Checking fort.90 ------------------------------" << std::endl;
		fort90fail = CheckFort90();
		std::cout << "---------------------------- End checking fort.90 ----------------------------" << std::endl;
		if(fort90fail)
		{
			std::cerr << "WARNING: fort.90 files do NOT match" << std::endl;
		}
		else
		{
			std::cout << "fort.90: PASS" << std::endl;
		}
	}

	//Look at STF
	if(STF)
	{
		std::cout << "------------------------ Checking singletrackfile.dat ------------------------" << std::endl;
		STFfail = CheckSTF();
		std::cout << "---------------------- End checking singletrackfile.dat ----------------------" << std::endl;
		if(STFfail)
		{
			std::cerr << "WARNING: singletrackfile.dat files do NOT match" << std::endl;
		}
		else
		{
			std::cout << "STF: PASS" << std::endl;
		}
	}
	//Look at extra_checks.txt
	ExtraChecksfail = PerformExtraChecks();

	std::cout << "---------------------------------- SUMMARY -----------------------------------" << std::endl;
	if(fort10)
	{
		if(fort10fail)
		{
			std::cout << "fort.10 DOES NOT MATCH" << std::endl;
		}
		else
		{
			std::cout << "fort.10 MATCHES" << std::endl;
		}
	}

	if(fort90)
	{
		if(fort90fail)
		{
			std::cout << "fort.90 DOES NOT MATCH" << std::endl;
		}
		else
		{
			std::cout << "fort.90 MATCHES" << std::endl;
		}
	}

	if(STF)
	{
		if(STFfail)
		{
			std::cout << "singletrackfile.dat DOES NOT MATCH" << std::endl;
		}
		else
		{
			std::cout << "singletrackfile.dat MATCHES" << std::endl;
		}
	}
	std::cout << "------------------------------------ EXIT ------------------------------------" << std::endl;

	//or together any fail bits.
	//If all tests pass this will return 0 (good)
	//if not we get something else out (bad)
	return (fort10fail || fort90fail || STFfail || ExtraChecksfail);
}

//Run the actual sixtrack binary
void RunSixTrack(char* argv[], int* Status)
{
	*Status = execl(argv[1], argv[1], (char*) 0);
}

/**
* Kills the running SixTrack process for Checkpoint/Resume (CR) builds
* @parm KillTime The time in seconds to wait until the run is killed
* @parm sixpid The pid of the process to kill
*/
void KillSixTrack(size_t KillTime, pid_t sixpid)
{
	//std::this_thread::sleep_for(std::chrono::seconds(KillTime));
	sleep(KillTime);
	std::cout << "KillSixTrack() - killing pid: " << sixpid << std::endl;
	int res = kill(sixpid, SIGKILL);
	std::cout << "KillSixTrack() kill result: " << res << std::endl;
}

/**
* Copies one input file to an output file.
* @parm InputFileName The name of the input file
* @parm OutputFileName The name of the output file
* @return A bool set to true if the operation succeeded
*/
bool CopyFile(std::string InputFileName, std::string OutputFileName)
{
	std::ifstream Input(InputFileName.c_str(), std::ios::binary);
	std::ofstream Output(OutputFileName.c_str(), std::ios::binary);

	if(!Input.good())
	{
		std::cerr << "ERROR: Could not open " << InputFileName << std::endl;
		return false;
	}

	if(!Input.good())
	{
		std::cerr << "ERROR: Could not open " << OutputFileName << std::endl;
		return false;
	}

	Output << Input.rdbuf();
	return true;
}

/**
* Checks if two fort.10 files are equivalent
* @return true if there is a problem, false if there is not
*/
bool CheckFort10()
{
	//return false if all is good, true if anything else happens

	//In this folder we should have fort.10 and a fort.10.canonical
	bool input1 = CopyFile("fort.10", "fort.20");
	bool input2 = CopyFile("fort.10.canonical", "fort.21");

	if(!input1 || !input2)
	{
		std::cerr << "WARNING: Could not perform fort.10 comparison" << std::endl;
		return true;
	}

	//Now again we fork() and exec()
	pid_t CheckFort10pid = fork();
	if(CheckFort10pid == -1)
	{
		std::cerr << "ERROR: Could not fork to start checkf10" << std::endl;
		return true;
	}

	//Check fork() status
	if(CheckFort10pid == 0)
	{
		//child, run checkf10
		int status = execl("./checkf10", "checkf10", (char*) 0);
		if(status == -1)
		{
			perror("ERROR: Could not execute checkf10");
			exit(EXIT_FAILURE);
		}
	}
	else
	{
		//main thread, wait()
		int waitpidStatus;
		waitpid(CheckFort10pid, &waitpidStatus, WUNTRACED);
		std::cout << "checkf10 finished running and returned: " << waitpidStatus << std::endl;
		return waitpidStatus;
	}
}

/**
* Checks if two fort.90 files are equivalent
* @return true if there is a problem, false if there is not
*/
bool CheckFort90()
{
	//return false if all is good, true if anything else happens
	//First we call read90 on each file.
	//Then we must do a binary comparison

	//Now again we fork() and exec()
	//Do this for the first file
	pid_t CheckFort90pid = fork();
	if(CheckFort90pid == -1)
	{
		std::cerr << "ERROR: Could not fork to start read90" << std::endl;
		return false;
	}

	//Check fork() status
	if(CheckFort90pid == 0)
	{
		//child, run read90
		int status = execl("./read90", "read90", "--fname", "fort.90", "--ofname", "fort.90.out", (char*) 0);
		if(status == -1)
		{
			perror("ERROR: Could not execute read90");
			exit(EXIT_FAILURE);
		}
	}
	else
	{
		//main thread, wait()
		int waitpidStatus;
		waitpid(CheckFort90pid, &waitpidStatus, WUNTRACED);
		std::cout << "read90 finished running on fort.90: " << waitpidStatus << std::endl;
		if(waitpidStatus != 0)
		{
			std::cerr << "ERROR: Problem running read90" << std::endl;
			return waitpidStatus;
		}
	}

	//Do this for the second file
	CheckFort90pid = fork();
	if(CheckFort90pid == -1)
	{
		std::cerr << "ERROR: Could not fork to start read90" << std::endl;
		return false;
	}

	//Check fork() status
	if(CheckFort90pid == 0)
	{
		//child, run read90
		int status = execl("./read90", "read90", "--fname", "fort.90.canonical", "--ofname", "fort.90.canonical.out", (char*) 0);
		if(status == -1)
		{
			perror("ERROR: Could not execute read90");
			exit(EXIT_FAILURE);
		}
	}
	else
	{
		//main thread, wait()
		int waitpidStatus;
		waitpid(CheckFort90pid, &waitpidStatus, WUNTRACED);
		std::cout << "read90 finished running on fort.90.canonical: " << waitpidStatus << std::endl;
		if(waitpidStatus != 0)
		{
			std::cerr << "ERROR: Problem running read90" << std::endl;
			return waitpidStatus;
		}
	}

	return !FileComparison("fort.90.out", "fort.90.canonical.out");

}

/**
* Checks if two singletrackfile.dat files are equivalent
* @return true if there is a problem, false if there is not
*/
bool CheckSTF()
{
	//return false if all is good, true if anything else happens
	//First we call read90 on each file.
	//Then we must do a binary comparison

	//Now again we fork() and exec()
	//Do this for the first file
	pid_t CheckSTFpid = fork();
	if(CheckSTFpid == -1)
	{
		std::cerr << "ERROR: Could not fork to start read90" << std::endl;
		return false;
	}

	//Check fork() status
	if(CheckSTFpid == 0)
	{
		//child, run read90
		int status = execl("./read90", "read90", "--STF", "--fname", "singletrackfile.dat", "--ofname", "singletrackfile.dat.out", (char*) 0);
		if(status == -1)
		{
			perror("ERROR: Could not execute read90");
			exit(EXIT_FAILURE);
		}
	}
	else
	{
		//main thread, wait()
		int waitpidStatus;
		waitpid(CheckSTFpid, &waitpidStatus, WUNTRACED);
		std::cout << "read90 finished running on singletrackfile.dat: " << waitpidStatus << std::endl;
		if(waitpidStatus != 0)
		{
			std::cerr << "ERROR: Problem running read90" << std::endl;
			return waitpidStatus;
		}
	}

	//Do this for the second file
	CheckSTFpid = fork();
	if(CheckSTFpid == -1)
	{
		std::cerr << "ERROR: Could not fork to start read90" << std::endl;
		return false;
	}

	//Check fork() status
	if(CheckSTFpid == 0)
	{
		//child, run read90
		int status = execl("./read90", "read90", "--STF", "--fname", "singletrackfile.dat.canonical", "--ofname", "singletrackfile.dat.canonical.out", (char*) 0);
		if(status == -1)
		{
			perror("ERROR: Could not execute read90");
			exit(EXIT_FAILURE);
		}
	}
	else
	{
		//main thread, wait()
		int waitpidStatus;
		waitpid(CheckSTFpid, &waitpidStatus, WUNTRACED);
		std::cout << "read90 finished running on singletrackfile.dat.canonical: " << waitpidStatus << std::endl;
		if(waitpidStatus != 0)
		{
			std::cerr << "ERROR: Problem running read90" << std::endl;
			return waitpidStatus;
		}
	}

	return !FileComparison("singletrackfile.dat.out", "singletrackfile.dat.canonical.out");

}
/**
* File comparison
*/
bool FileComparison(std::string FileName1, std::string FileName2)
{
	std::ifstream f1(FileName1.c_str(), std::ios::binary);
	std::ifstream f2(FileName2.c_str(), std::ios::binary);

	size_t length1 = 0;
	size_t length2 = 0;
	if(f1.good())
	{
		f1.seekg(0,f1.end);
		length1 = f1.tellg();
	}
	else
	{
		std::cout << "Could not open " << FileName2 << std::endl;
		return false;
	}

	if(f2.good())
	{
		f2.seekg(0,f2.end);
		length2 = f2.tellg();
	}
	else
	{
		std::cout << "Could not open " << FileName2 << std::endl;
		return false;
	}

	std::cout << "size " << FileName1 << ": " << length1 << std::endl;
	std::cout << "size " << FileName2 << ": " << length2 << std::endl;

	if(length1 != length2)
	{
		return false;
	}

	//Same size files, now to check the contents
	const size_t BufferSize = 1024;
	char f1Buffer[BufferSize];
	char f2Buffer[BufferSize];

	f1.seekg(0);
	f2.seekg(0);

	size_t position = 0;

	while(position <= length1)
	{
		size_t rsize = std::min(length1, BufferSize);

		f1.read(f1Buffer,rsize);
		f2.read(f2Buffer,rsize);

		//std::cerr << FileName1 << " is different from " << FileName2 << " at " << length2-length1 << " - comparison: " << comparison << std::endl;

		int comparison = memcmp(f1Buffer, f2Buffer, rsize);
		if(comparison != 0)
		{
			std::cerr << FileName1 << " is different from " << FileName2 << " at " << position << " - memcmp() comparison: " << comparison << std::endl;
			return false;
		}

		position+= rsize;
	}
	return true;
}

bool PerformExtraChecks()
{
	std::cout << "--------------------------- Performing extra checks ---------------------------" << std::endl;
	bool AllTests = false;
	std::ifstream extra_checks_in("extra_checks.txt");
	if(extra_checks_in.good())
	{

		std::cout << "Opened extra_checks.txt" << std::endl;
		//Format should be some file to check followed by a command
		while(extra_checks_in.good())
		{
			std::string StringBuffer;
			std::string FileName;

			extra_checks_in >> FileName;
			getline(extra_checks_in, StringBuffer);
			if(FileName != "")
			{
				std::cout << "Performing extra checks on " << FileName << std::endl;
				bool ThisTest = !FileComparison(FileName, FileName + ".canonical");
				if(ThisTest)
				{
					std::cerr << "WARNING: Extra check on " << FileName << " failed!" << std::endl;
					AllTests = true;
				}
				else
				{
					std::cout << "Extra check on " << FileName << " MATCHES" << std::endl;
				}
			}
		}
	}
	else
	{
		std::cout << "Could not open extra_checks.txt" << std::endl;
	}

	std::cout << "------------------------------- End extra checks ------------------------------" << std::endl;
	return AllTests;
}

std::vector<int> ParseKillTimes(char* in)
{
	std::vector<int> KillTimes;
	std::string input = in;
	size_t pos = 0;
	size_t oldpos = 0;
	while(pos != std::string::npos)
	{
		pos = input.find(",", oldpos);
		int number = atoi(input.substr(oldpos,pos-oldpos).c_str());
		oldpos=pos+1;

		KillTimes.push_back(number);
	}

	std::cout << "Will try and kill CR run " << KillTimes.size() << " times." << std::endl;
	std::cout << "Killing after ";
	for(int k = 0; k < KillTimes.size(); k++)
	{
		std::cout << KillTimes.at(k) << "\t";
	}
	std::cout << " seconds." << std::endl;
	return KillTimes;
}
