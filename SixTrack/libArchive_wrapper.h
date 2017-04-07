#ifndef LIBARCHIVE_WRAPPER_H
#define LIBARCHIVE_WRAPPER_H

#ifdef __cplusplus
extern "C" {
#endif

void write_archive(const char* const outname, char** filename, int nFiles);
void list_archive(const char* const infile);
void list_archive_get(const char* const infile, char** filenames, int* nfiles, const int buffsize);
void read_archive(const char* const infile, const char* const extractFolder);

#ifdef __cplusplus
}
#endif

#endif
