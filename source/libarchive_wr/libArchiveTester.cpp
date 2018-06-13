//Little program to open and list archives.

#include "libArchive_wrapper.h"

#include <stdlib.h>
#include <iostream>
using namespace std;

int main(int argc, char* argv[]){
  if (argc != 2) {
    cout << "Usage: ./libArchiveTester path-to-zipfile" << endl;
    cout << "It will then print a list of the file contents." << endl;
    cout << "Afterwards it tries to unzip the archive; this may fail if the ZLIB archive was not correctly linked." << endl;
    return EXIT_FAILURE;
  }

  cout << "Listing the archive:" << endl;
  list_archive(argv[1]);

  cout << endl;

  cout << "OK, now inflating into folder 'libarchive_test' in current directory!" << endl;
  read_archive(argv[1], "libarchive_test");

  cout << "Done!" << endl;
  return EXIT_SUCCESS;
}
