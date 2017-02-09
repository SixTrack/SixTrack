#ifndef LIBARCHIVE_WRAPPER_H
#define LIBARCHIVE_WRAPPER_H

void write_archive(const char* const outname, char** filename, int nFiles);
void list_archive(const char* const infile);
void read_archive(const char* const infile, const char* extractFolder);

#endif
