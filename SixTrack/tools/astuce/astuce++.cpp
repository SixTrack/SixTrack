#include "astuce++.h"
#include <cctype>
#include <iomanip>
#include <cstdlib>

void Astuce::Run()
{
	LoadInputMask();
	ReadInputMaskFile();
	if(OutputFlag){PrintInputSettings();}
	ReadInputFiles();
	ParseSourceFile();
	DoOutput();
}

void Astuce::LoadInputMask()
{
	std::ifstream *infile = new std::ifstream(InputFileName);

	if(infile->good())
	{
		if(OutputFlag)
		{
			*LogStream << "Opened " << InputFileName << std::endl;
		}
	}
	else
	{
		std::cerr << "Could not open " << InputFileName << std::endl;
		exit(EXIT_FAILURE);
	}

	std::string c_buffer;
	while(infile->good())
	{
		getline(*infile, c_buffer);
		if(!infile->eof())
		{
			InputBuffer.push_back(c_buffer);
		}
	}

	infile->close();
	delete infile;
}

void Astuce::ReadInputMaskFile()
{
	size_t line = 1;
	std::string Buffer;

	std::vector<std::string>::const_iterator itr = InputBuffer.begin();
	while(itr != InputBuffer.end())
	{
		if(line == 1)
		{
			Buffer = *itr;

			size_t space = 0;
			size_t old_pos = 0;
			while(space!=std::string::npos)
			{
				space = Buffer.find(" ", space+1);
				InputFileNames.insert(Buffer.substr(old_pos,space-old_pos));
				old_pos = space+1;
			}
		}
		else if(line == 2)
		{
			Buffer = *itr;

			size_t space = Buffer.find_first_of(" \n", 0);
			OutputFileName = Buffer.substr(0,space);

		}
		else
		{
			if(itr->substr(0,2) == "EX" || itr->substr(0,2) == "ex")
			{
				//Exit processing
				break;
			}
			else if(itr->substr(0,2) == "DF" || itr->substr(0,2) == "df")
			{
				//flag defines
				//Need to extract first, then split on a ',' char
				size_t start = itr->find_first_not_of(" ", 2);
				if(start!=std::string::npos)
				{
					std::string defines = itr->substr(start,itr->size()-start);
					size_t step = 0;
					size_t old_step = 0;
					while(step!=std::string::npos)
					{
						step = defines.find(",", step+1);
						Flags.insert(FixCase(defines.substr(old_step,step-old_step)));
						old_step = step+1;
					}
				}
			}
			else if(itr->substr(0,2) == "E " || itr->substr(0,2) == "e ")
			{
				//Need to extract first, then split on a ',' char
				size_t start = itr->find_first_not_of(" ", 1);
				if(start!=std::string::npos)
				{
					std::string defines = itr->substr(start,itr->size()-start);
					size_t step = 0;
					size_t old_step = 0;
					while(step!=std::string::npos)
					{
						step = defines.find(",", step+1);
						Decks.insert(FixCase(defines.substr(old_step,step-old_step)));
						old_step = step+1;
					}
				}
			}
			else
			{
				std::cerr << "ERROR processing " << InputFileName << std::endl;
				std::cerr << "ERROR: Bad entry in input: " << *itr << std::endl;
				exit(EXIT_FAILURE);
			}
		}

		line++;
		itr++;
	}
}

std::string Astuce::FixCase(std::string in)
{
	std::string out;
	for(size_t n=0; n < in.size(); n++)
	{
		out.push_back(toupper(in.at(n)));
	}
	return out;
}

void Astuce::PrintInputSettings()
{
	*LogStream << "\nINPUT FILES:" << "\n";
	std::set<std::string>::const_iterator InputFileNames_itr = InputFileNames.begin();
	while(InputFileNames_itr != InputFileNames.end())
	{
		*LogStream  << *InputFileNames_itr << "\n";
		InputFileNames_itr++;
	}

	*LogStream << "\nOUTPUT FILE:\n" << OutputFileName << "\n";

	*LogStream << "\nSELECTED FLAGS:" << "\n";
	std::set<std::string>::const_iterator Flags_itr = Flags.begin();
	while(Flags_itr != Flags.end())
	{
		*LogStream << *Flags_itr << "\n";
		Flags_itr++;
	}

	*LogStream << "\nSELECTED DECKS:" << std::endl;
	std::set<std::string>::const_iterator Decks_itr = Decks.begin();
	while(Decks_itr != Decks.end())
	{
		*LogStream << *Decks_itr << "\n";
		Decks_itr++;
	}
	*LogStream << std::endl;
}

void Astuce::ReadInputFiles()
{
	std::set<std::string>::const_iterator InputFileNames_itr =  InputFileNames.begin();
	while(InputFileNames_itr != InputFileNames.end())
	{
		std::ifstream* InputFileRead = new std::ifstream(InputFileNames_itr->c_str(), std::ios::binary);
		if (!InputFileRead->good())
		{
			std::cerr << "ERROR processing " << InputFileName << std::endl;
			std::cerr << "ERROR: Error reading input file '" << InputFileNames_itr->c_str() << "'"<< std::endl;
			exit(EXIT_FAILURE);
		}
		char cbuffer;
		std::streambuf* sbuffer = InputFileRead->rdbuf();

		//Make sure every file in the buffer begins with a \n,
		// which makes looking for +cd,+dk etc. easier.
		// It is also helpfull when reading multiple source files
		// in case one of them is ending on a newline.
		InputFileBufferStream << '\n';
		
		while(sbuffer->sgetc() != EOF)
		{
			if((cbuffer = sbuffer->sbumpc()))
			{
				InputFileBufferStream << cbuffer;
			}
		}

		InputFileRead->close();
		delete InputFileRead;

		InputFileNames_itr++;
	}
	//*LogStream << InputFileBufferStream.str();
	InputFileBuffer = InputFileBufferStream.str();
}

void Astuce::ParseSourceFile()
{
//Language
//+cd "COMDECK"
//+ca "CALL"
//+if if
//.NOT. not
//.AND. and
//.OR. or
//+ei end if
//+el else if (not used in sixtrack?)
//+dk "DECK"

	//Step 1: Extract all +cd blocks into a string map
	ExtractCalls();

	//Step 2: Extract all +dk blocks into a string map
	ExtractDecks();

	//Step 3: Verify the input - do the requested decks exist?
	VerifyInput();

	//Step 4: enable/disable +if/+ei blocks
	ExpandIF();

	//Step 5: Insert +cd blocks into the code where +ca exists
	ExpandCA();
}

void Astuce::ExtractCalls()
{
	size_t CurrentPosition = 0;
	size_t LineEnd = 0;
	size_t Next_cd = 0;
	size_t Next_dk = 0;

	while(CurrentPosition < InputFileBuffer.size())
	{
		CurrentPosition = InputFileBuffer.find("\n+cd", CurrentPosition);

		if(CurrentPosition == std::string::npos)
		{
			break;
		}

		LineEnd = InputFileBuffer.find("\n", CurrentPosition+1);
		size_t NameStart = InputFileBuffer.find_first_not_of(" ", CurrentPosition+4);
		size_t NameEnd = InputFileBuffer.find_first_of(" !\n", NameStart+1);
		std::string BlockName = FixCase(InputFileBuffer.substr(NameStart,NameEnd-NameStart));

		//This could end at either the next +cd or the next +dk
		//Find where the next one of these takes place and read until there.
		Next_cd = InputFileBuffer.find("\n+cd", LineEnd);
		Next_dk = InputFileBuffer.find("\n+dk", LineEnd);
		CurrentPosition = std::min(Next_cd, Next_dk);

		std::string BlockEntry = InputFileBuffer.substr(LineEnd+1, CurrentPosition-LineEnd-1);

		//Now split this +cd block into lines
		std::list<LineStorage> Lines;
		LineStorage Entry;

		size_t StartPosition = 0;
		size_t BlockPosition = 0;

		while(BlockPosition < BlockEntry.size())
		{
			BlockPosition = BlockEntry.find("\n", StartPosition);
			if(BlockPosition == std::string::npos)
			{
				BlockPosition = BlockEntry.size();
			}

			Entry.Text = BlockEntry.substr(StartPosition, BlockPosition-StartPosition);
			Entry.Enabled = true;

			Lines.push_back(Entry);
			StartPosition = BlockPosition+1;
		}

		std::pair<std::map<std::string, std::list<LineStorage> >::const_iterator,bool > check;
		check = CallStorage.insert(std::make_pair(BlockName,Lines));

		if(check.second == false)
		{
			std::cerr << "ERROR processing " << InputFileName << std::endl;
			std::cerr << "Multiple definitions of call: \"" << BlockName << "\"" << std::endl;
			exit(EXIT_FAILURE);
		}
	}

	if(OutputFlag)
	{
		std::map<std::string, std::list<LineStorage> >::const_iterator it = CallStorage.begin();
		while(it!=CallStorage.end())
		{
			*LogStream << it->first << " ";
			it++;
		}
		*LogStream << std::endl;
		*LogStream << "Total \"calls\": " << CallStorage.size() << "\n" << std::endl;
	}
}

void Astuce::ExtractDecks()
{
	size_t CurrentPosition = 0;
	size_t LineEnd = 0;
	size_t Next_cd = 0;
	size_t Next_dk = 0;

	while(CurrentPosition < InputFileBuffer.size())
	{
		CurrentPosition = InputFileBuffer.find("\n+dk", CurrentPosition);
		
		if(CurrentPosition == std::string::npos)
		{
			break;
		}

		LineEnd = InputFileBuffer.find("\n", CurrentPosition+1);
		size_t NameStart = InputFileBuffer.find_first_not_of(" ", CurrentPosition+4);
		size_t NameEnd = InputFileBuffer.find_first_of(" !\n", NameStart+1);
		std::string BlockName = FixCase(InputFileBuffer.substr(NameStart,NameEnd-NameStart));

		//This could end at either the next +cd or the next +dk
		//Find where the next one of these takes place and read until there.
		Next_cd = InputFileBuffer.find("\n+cd", LineEnd);
		Next_dk = InputFileBuffer.find("\n+dk", LineEnd);
		CurrentPosition = std::min(Next_cd, Next_dk);

		std::string BlockEntry = InputFileBuffer.substr(LineEnd+1, CurrentPosition-LineEnd-1);

		//Now split this +dk block into lines
		std::list<LineStorage> Lines;
		LineStorage Entry;

		size_t StartPosition = 0;
		size_t BlockPosition = 0;

		while(BlockPosition < BlockEntry.size())
		{
			BlockPosition = BlockEntry.find("\n", StartPosition);
			if(BlockPosition == std::string::npos)
			{
				BlockPosition = BlockEntry.size();
			}

			Entry.Text = BlockEntry.substr(StartPosition, BlockPosition-StartPosition);
			Entry.Enabled = true;

			Lines.push_back(Entry);
			StartPosition = BlockPosition+1;
		}

		std::pair<std::map<std::string, std::list<LineStorage> >::const_iterator,bool > check;
		check = DeckStorage.insert(std::make_pair(BlockName,Lines));

		if(check.second == false)
		{
			std::cerr << "ERROR processing " << InputFileName << std::endl;
			std::cerr << "Multiple definitions of deck: \"" << BlockName << "\"" << std::endl;
			exit(EXIT_FAILURE);
		}
	}

	if(OutputFlag)
	{
		std::map<std::string, std::list<LineStorage> >::const_iterator it = DeckStorage.begin();
		while(it!=DeckStorage.end())
		{
			*LogStream << it->first << " ";
			it++;
		}
		*LogStream << std::endl;
		*LogStream << "Total \"decks\": " << DeckStorage.size() << "\n" << std::endl;
	}
}

void Astuce::ExpandIF()
{
	//start with +cd
	//Remember we can have nested if levels!

	std::map<std::string, std::list<LineStorage> >::iterator CallStorage_itr = CallStorage.begin();
	//This loops over each +cd block
	while(CallStorage_itr != CallStorage.end())
	{
		DoExpandIF(CallStorage_itr);
		CallStorage_itr++;
	}

	//Then +dk
	std::map<std::string, std::list<LineStorage> >::iterator DeckStorage_itr = DeckStorage.begin();
	//This loops over each +dk block
	while(DeckStorage_itr != DeckStorage.end())
	{
		DoExpandIF(DeckStorage_itr);
		DeckStorage_itr++;
	}
}

void Astuce::DoExpandIF(std::map<std::string, std::list<LineStorage> >::iterator Storage_itr)
{
	std::string BlockName = Storage_itr->first;

	if(OutputFlag)
	{
		*LogStream << "Processing: " << BlockName << " with " << Storage_itr->second.size() << " lines." << std::endl;
	}

	//This keeps track of the current state of +if/+ei statements, and if they are enabled or not.
	EnabledStates.push(true);

	std::list<LineStorage>::iterator Lines_itr = Storage_itr->second.begin();
	while(Lines_itr != Storage_itr->second.end())
	{
		ProcessIfBlock(Lines_itr);
		Lines_itr++;
	}

	//Check the current state
	if(EnabledStates.size() != 1)
	{
		std::cerr << "ERROR processing " << InputFileName << std::endl;
		std::cerr << "ERROR: Mismatched number of +if and +ei statements, check your .s files!" << std::endl;
		std::cerr << "Currently processing section " << BlockName << std::endl;
		std::cerr << "Current level of +if/+ei: " << EnabledStates.size() << " (Expected: 1)" << std::endl;
		exit(EXIT_FAILURE);
	}
	EnabledStates.pop();


	if(OutputFlag)
	{
		Lines_itr = Storage_itr->second.begin();
		while(Lines_itr != Storage_itr->second.end())
		{
			*LogStream << BlockName << "\t" << Lines_itr->Enabled << " " << Lines_itr->Text << std::endl;
			Lines_itr++;
		}

	}
}

void Astuce::ProcessIfBlock(std::list<LineStorage>::iterator line)
{
	//Re-evaluate the enabled state after each +if and +ei statement
	//Before each +if block we need to store the previous state and restore at +ei (returning to the same depth)
	//If a parent block is disabled, all children should also be disabled.
	//Done with a std::stack<bool> and push()/pop()/top()

	//Does this block have a +if statement?
	if(line->Text.substr(0,3) == "+if")
	{
		if(line->Text.substr(3,1) != " ")
		{
			std::cerr << "ERROR processing \"" << InputFileName << "\"" << std::endl;
			std::cerr << "+if statement with no space following \"+if\" : \"" << line->Text << "\"" << std::endl;
			exit(EXIT_FAILURE);
		}
		//Disable this line
		line->Enabled = false;

		//The result of the +if
		bool EnabledState = false;

		//Find the name of this statement
		size_t NameStart = line->Text.find_first_not_of(" ", 4);
		if(NameStart == std::string::npos)
		{
			std::cerr << "ERROR processing \"" << InputFileName << "\"" << std::endl;
			std::cerr << "+if statement with no associated flag? : \"" << line->Text << "\"" << std::endl;
			exit(EXIT_FAILURE);
		}

		//Strip off comments
		size_t NameEnd = line->Text.find_first_of("!", NameStart+1);
		//Strip off trailing whitespace
		if (NameEnd != std::string::npos) NameEnd-=1;
		NameEnd = line->Text.find_last_not_of(" ",NameEnd);
		if (NameEnd != std::string::npos) NameEnd+=1;
		
		std::string FlagName = FixCase(line->Text.substr(NameStart,NameEnd-NameStart));

		if (FlagName.find(" ") != std::string::npos)
		{
			std::cerr << "ERROR processing \"" << InputFileName << "\"" << std::endl;
			std::cerr << "The +if statement '" << line->Text << "' contains a space inside the actual statement '" << FlagName << "'." << std::endl;
			exit(EXIT_FAILURE);
		}
		
		//Check that there are no spaces inside the STATEMENT
		if(OutputFlag)
		{
			*LogStream << "Evaluating flag: '" << FlagName << "'" << std::endl;
		}

		
		size_t Position = 0;
		size_t NextOR = 0;
		size_t NextAND = 0;
		std::vector<size_t> Splits;
		std::vector<std::string> FlagNames;
		std::vector<bool> And;
		std::vector<bool> Or;
		std::vector<bool> States;

		Splits.push_back(Position);
		while(NextOR != std::string::npos && NextAND != std::string::npos)
		{
			NextOR = FlagName.find(".OR.", Position);
			NextAND = FlagName.find(".AND.", Position);
			Position = std::min(NextOR, NextAND);
			if(Position != std::string::npos)
			{
				Splits.push_back(Position);
			}
			Position++;
		}
		Splits.push_back(FlagName.size());

		for(size_t n=1; n < Splits.size(); n++)
		{
			std::string SplitSubString = FlagName.substr(Splits.at(n-1),Splits.at(n));
			FlagNames.push_back(SplitSubString);
		}

		for(size_t n=0; n < FlagNames.size(); n++)
		{
			bool NOT = false;
			bool AND = false;
			bool OR = false;

			if(FlagNames.at(n).substr(0,5) == ".AND.")
			{
				AND = true;
				FlagNames.at(n).erase(0,5);
			}
			else if(FlagNames.at(n).substr(0,4) == ".OR.")
			{
				OR = true;
				FlagNames.at(n).erase(0,4);
			}
			if(FlagNames.at(n).substr(0,5) == ".NOT.")
			{
				FlagNames.at(n).erase(0,5);
				NOT = true;
			}

			//Check that there is actually a flag following the .AND./.OR./.NOT.
			if (FlagNames.at(n).length() == 0) {
				std::cerr << "ERROR processing \"" << InputFileName << "\"" << std::endl;
				std::cerr << "The +if statement '" << line->Text << "' has an orphan .AND./.OR./.NOT." << std::endl;
				exit(EXIT_FAILURE);
			}
			
			//Enable/disable lines
			std::set<std::string>::const_iterator Flags_itr = Flags.find(FlagNames.at(n));
			if(Flags_itr != Flags.end())
			{
				//flag is found -> enabled
				if(NOT == false){EnabledState = true;}
				else{EnabledState = false;}
			}
			else
			{
				//flag is NOT found -> enabled
				//.not.xxx and xxx is not enabled
				if(NOT == false){EnabledState = false;}
				else{EnabledState = true;}
			}
			Or.push_back(OR);
			And.push_back(AND);
			States.push_back(EnabledState);
		}

		//Do the first
		EnabledState = States.at(0);

		for(size_t n = 1; n < States.size(); n++)
		{
			if(And.at(n) == true)
			{
				EnabledState = (EnabledState & States.at(n));
			}
			else if(Or.at(n) == true)
			{
				EnabledState = (EnabledState | States.at(n));
			}
			else
			{
				//Should never hit this
				std::cerr << "ERROR processing \"" << InputFileName << "\"" << std::endl;
				std::cerr << "ERROR: in .if. .and. .or. processing - Multiple entries and neither .and. or .or. true at entry " << n << std::endl;
				exit(EXIT_FAILURE);
			}
		}

		if(EnabledStates.empty())
		{
			std::cerr << "ERROR processing \"" << InputFileName << "\"" << std::endl;
			std::cerr << "ERROR: Mismatched number of +if and +ei statements, check your .s files!" << std::endl;
			std::cerr << "Currently processing \"" << line->Text << "\"" << std::endl;
			exit(EXIT_FAILURE);
		}
		//AND (&) the result of the above with the current stack top
		EnabledStates.push(EnabledState & EnabledStates.top());
	}
	else if(line->Text.substr(0,3) == "+ei")
	{
		//Disable this line
		line->Enabled = false;

		if(EnabledStates.empty())
		{
			std::cerr << "ERROR processing \"" << InputFileName << "\"" << std::endl;
			std::cerr << "ERROR: Mismatched number of +if and +ei statements, check your .s files!" << std::endl;
			std::cerr << "Currently processing " << line->Text << std::endl;
			exit(EXIT_FAILURE);
		}
		EnabledStates.pop();
	}
	else
	{
		//If a general line is to be enabled or disabled, this depends on the current state
		if(EnabledStates.empty())
		{
			std::cerr << "ERROR processing \"" << InputFileName << "\"" << std::endl;
			std::cerr << "ERROR: Mismatched number of +if and +ei statements, check your .s files!" << std::endl;
			std::cerr << "Currently processing " << line->Text << std::endl;
			exit(EXIT_FAILURE);
		}
		line->Enabled = EnabledStates.top();
		return;
	}
}

void Astuce::ExpandCA()
{
	bool Replaced = true;

	//First step, go through existing +cd blocks and check for +ca
	while(Replaced)
	{
		//Repeat until no +ca exists
		Replaced = false;

		std::map<std::string, std::list<LineStorage> >::iterator CallStorage_itr = CallStorage.begin();

		//This loops over each +cd block
		while(CallStorage_itr != CallStorage.end())
		{
			//Perform replace
			if(ReplaceCA(CallStorage_itr))
			{
				Replaced = true;
			}

			CallStorage_itr++;
		}
	}

	//Add +ca blocks where needed into +dk
	std::map<std::string, std::list<LineStorage> >::iterator DeckStorage_itr = DeckStorage.begin();
	//This loops over each +cd block
	while(DeckStorage_itr != DeckStorage.end())
	{
		ReplaceCA(DeckStorage_itr);
		DeckStorage_itr++;
	}
}

bool Astuce::ReplaceCA(std::map<std::string, std::list<LineStorage> >::iterator Storage_itr)
{
	//Storage_itr = Each block
	//Lines_itr = Each line in the block

	//CallStorage_itr A selected global block to insert
	//Insert_itr Each line in the global block

	bool Replaced = false;

	std::string BlockName = Storage_itr->first;
	
	std::list<LineStorage>::iterator Lines_itr = Storage_itr->second.begin();
	while(Lines_itr != Storage_itr->second.end())
	{
		if(Lines_itr->Enabled == true)
		{
			if(Lines_itr->Text.substr(0,3) == "+ca")
			{
				//Get the +ca name
				size_t NameStart = Lines_itr->Text.find_first_not_of(" ", 4);
				size_t NameEnd = Lines_itr->Text.find_first_of(" !", NameStart+1);
				std::string CallName = FixCase(Lines_itr->Text.substr(NameStart,NameEnd-NameStart));

				//Find this name
				std::map<std::string, std::list<LineStorage> >::iterator CallStorage_itr = CallStorage.begin();
				CallStorage_itr = CallStorage.find(CallName);
				if(CallStorage_itr != CallStorage.end())
				{
					//Delete the +ca statement line
					Lines_itr = Storage_itr->second.erase(Lines_itr);

					LineStorage Title;
					Title.Enabled=true;
					Title.Text="!!! START " + CallName;
///					Storage_itr->second.insert(Lines_itr,Title);

					std::list<LineStorage>::iterator Insert_itr = CallStorage_itr->second.begin();
					while(Insert_itr != CallStorage_itr->second.end())
					{
						//Do the insert
						Storage_itr->second.insert(Lines_itr, *Insert_itr);

						//Advance to the next line
						Insert_itr++;
					}
					Replaced = true;


					Title.Text="!!! END " + CallName;
///					Storage_itr->second.insert(Lines_itr,Title);
				}
				else
				{
					std::cerr << "ERROR processing \"" << InputFileName << "\"" << std::endl;
					std::cerr << "Problematic line: \"" << Lines_itr->Text << "\"" << std::endl;
					std::cerr << "Searched for \"+cd " << CallName << "\" and it was not defined anywhere!" << std::endl;
					std::cerr << "The following +cd blocks exist:" << std::endl;
					CallStorage_itr = CallStorage.begin();
					while(CallStorage_itr != CallStorage.end())
					{
						std::cerr << "\"" << CallStorage_itr->first << "\" ";
						CallStorage_itr++;
					}
					std::cerr << std::endl;

					std::cerr << "Calls entry count: " << CallStorage.size() << std::endl;
					std::cerr << "Decks entry count: " << DeckStorage.size() << std::endl;

					exit(EXIT_FAILURE);
				}

				Lines_itr--;
			}
		}
		Lines_itr++;
	}
	return Replaced;
}

void Astuce::DoOutput()
{
	std::set<std::string>::const_reverse_iterator Decks_itr = Decks.rbegin();
	while(Decks_itr != Decks.rend())
	{
		std::map<std::string, std::list<LineStorage> >::iterator DeckStorage_itr = DeckStorage.begin();

		DeckStorage_itr = DeckStorage.find(*Decks_itr);
		if(DeckStorage_itr != DeckStorage.end())
		{
			std::list<LineStorage>::iterator Lines_itr = DeckStorage_itr->second.begin();
			while(Lines_itr != DeckStorage_itr->second.end())
			{
				if(Lines_itr->Enabled == true)
				{
//					OutputFileBufferStream << Lines_itr->Text.substr(0,80) << "\n";
					OutputFileBufferStream << Lines_itr->Text << "\n";
				}
				Lines_itr++;
			}
		}

		Decks_itr++;
	}

	std::ofstream* OutputStream = new std::ofstream(OutputFileName.c_str());
	*OutputStream << OutputFileBufferStream.rdbuf();
	OutputStream->close();
	delete OutputStream;
}

void Astuce::SetOutput(bool flg)
{
	OutputFlag = flg;
}

void Astuce::LogToFile(bool flg)
{
	if(OutputFlag)
	{
		if(flg)
		{
			std::string fname = (std::string)InputFileName + ".log";
			LogStream = new std::ofstream(fname.c_str());
		}
	}
}

void Astuce::VerifyInput()
{
	std::set<std::string>::const_iterator Decks_itr = Decks.begin();
	while(Decks_itr != Decks.end())
	{
		std::map<std::string, std::list<LineStorage> >::iterator DeckStorage_itr = DeckStorage.begin();
		DeckStorage_itr = DeckStorage.find(*Decks_itr);
		if(DeckStorage_itr == DeckStorage.end())
		{
			std::cerr << "ERROR processing \"" << InputFileName << "\"" << std::endl;
			std::cerr << "ERROR: Deck \"" << *Decks_itr << "\" requested but was not found in the input files!" << std::endl;
			exit(EXIT_FAILURE);
		}

		Decks_itr++;
	}
}
