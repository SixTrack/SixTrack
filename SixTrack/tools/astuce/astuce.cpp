#include <iostream>
#include "astuce++.h"

/**
* Line 1: A list of input file names
* Line 2: The output file name
* Line 3: df (flags)
* line 4+: extract decks
* line n: ex (exit)
*/
int main(int argc, char* argv[])
{
	if(argc != 2)
	{
		std::cout << "Usage: ./astuce <input mask name>" << std::endl;
		return EXIT_FAILURE;
	}

	Astuce* ast = new Astuce(argv[1]);
	ast->SetOutput(true);
	ast->LogToFile(true);
	ast->Run();
	delete ast;

	return EXIT_SUCCESS;
}
