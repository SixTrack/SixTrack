#include <iostream>
#include <fstream>

#include <vector>
#include <stack>
#include <set>
#include <map>
#include <list>

#include <string>
#include <sstream>

class Astuce
{
public:

	Astuce(char* in, bool LogFlag=false) : InputFileName(in), OutputFlag(LogFlag), LogStream(&std::cout) {}

	~Astuce()
	{
		if(LogStream != &std::cout)
		{
			dynamic_cast<std::ofstream*>(LogStream)->close();
			delete LogStream;
		}
	}

	class LineStorage
	{
		public:
		bool Enabled;
		std::string Text;
	};

	//The input ast file
	char* InputFileName;
	std::vector<std::string> InputBuffer;

	//The input file names to load
	std::set<std::string> InputFileNames;
	std::set<std::string> Flags;
	std::set<std::string> Decks;
	std::string OutputFileName;

	std::stringstream InputFileBufferStream;
	std::string InputFileBuffer;

	std::map<std::string, std::list<LineStorage> > CallStorage;
	std::map<std::string, std::list<LineStorage> > DeckStorage;

	//Current Enabled/Disabled state of code lines
	std::stack<bool> EnabledStates;

	std::stringstream OutputFileBufferStream;


	//Enable/disable stdout
	bool OutputFlag;

	//The stream to send the logging output to
	std::ostream* LogStream;

	//Sets the above Output flag
	void SetOutput(bool);

	//Logs output for each run to "<InputFileName>.log"
	//Normal default is std::cout
	void LogToFile(bool);


	//Do everything!
	void Run();

	//Loads up the input .ast file
	void LoadInputMask();

	//Parses the input .ast file
	void ReadInputMaskFile();

	//Print the input settings read from the input .ast file
	void PrintInputSettings();

	//Read the input .s files into a buffer
	void ReadInputFiles();

	//Does the first parse
	void ParseSourceFile();

	//Extracts +cd flags
	void ExtractCalls();

	//Extracts +dk decks
	void ExtractDecks();

	//Expands +if/+ei statements (enables/disables lines)
	void ExpandIF();

	//if processing for each "block"
	void DoExpandIF(std::map<std::string, std::list<LineStorage> >::iterator);

	//Sub-block processing
	void ProcessIfBlock(std::list<LineStorage>::iterator);

	//replace/insert +ca statements
	void ExpandCA();

	//replace/insert +ca statements
	bool ReplaceCA(std::map<std::string, std::list<LineStorage> >::iterator Storage_itr);

	//Write the output file
	void DoOutput();

	/* Helper functions */
	//Input is mixed upper/lower case string, fix this
	std::string FixCase(std::string);

};

