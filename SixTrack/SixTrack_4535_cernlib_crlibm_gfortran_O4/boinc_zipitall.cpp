#include <stdio.h>
#include <unistd.h>
#include <vector>
#include <string>
#include <cstdlib>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <time.h>
#include <stdlib.h>

#include "error_numbers.h"
#include "filesys.h"

#include "config.h"
#include "util.h"
#include "boinc_api.h"
#include "boinc_zip.h"

using std::vector;
using std::string;

extern "C" {

    void boinc_zipitall_()
    {
        ZipFileList plist;
        string zipfile = "Sixout.zip";

        bool
status=boinc_filelist(".","fort.",&plist,SORT_NAME|SORT_DESCENDING,true);
       
status=boinc_filelist(".","stderr.",&plist,SORT_NAME|SORT_DESCENDING,false);
        boinc_zip(ZIP_IT,zipfile,&plist);
        return;
    }

   void boinc_unzip_()
   {
        string filepath;
        int retval = boinc_resolve_filename_s("Sixin.zip",filepath);
        boinc_zip(UNZIP_IT,filepath,NULL);
   }

}
