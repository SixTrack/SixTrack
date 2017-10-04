//Little program to open and list archives.

#include "libArchive_wrapper.h"

#include <stdlib.h>
#include <iostream>
using namespace std;

int main(int argc, char* argv[]){
  if (argc != 2) {
    cout << "Usage: ./libArchiveTester path-to-zipfile" << endl;
    cout << "It will then print a list of the file contents." << endl;
    return EXIT_FAILURE;
  }

  list_archive(argv[1]);
  cout << "OK, now inflating!"<<endl;
  read_archive(argv[1], "libarchive_test");

  return EXIT_SUCCESS;
}
