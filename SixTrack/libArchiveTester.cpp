//Little program to open and list archives.

#include "libArchive_wrapper.h"

#include <iostream>
using namespace std;

int main(int argc, char* argv[]){
  if (argc != 2) {
    cout << "Usage: ./libArchiveTester path-to-zipfile" << endl;
    cout << "It will then print a list of the file contents." << endl;
    return 0;
  }

  list_archive(argv[1]);

  return 0;
}
