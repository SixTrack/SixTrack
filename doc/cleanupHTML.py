#!/usr/bin/env python3

#
# This script cleans up the html file(s) necessary for the /SixTrack/web/docs folder on the website.
# The files are written again as php files with the necessary headers and includes for the website.
# Written by Veronica Berglyd Olsen, Feb 2018
#

import sys

from os import path, listdir

def wrapHTML(htmlDir, htmlFile, fileList):
    
    phpFile  = htmlFile[0:-5]+".php"
    htmlBody = ""
    htmlHead = ""
    genStr   = ""
    whereNow = 0
    
    print("")
    print("Loading HTML File ...")
    
    with open(path.join(htmlDir,htmlFile),encoding="utf-8") as inFile:
        for theLine in inFile:
            
            theLine = theLine.rstrip()
            if theLine == "": continue
            
            if theLine[0:6] == "<head>":
                whereNow = 1
                continue
            
            if theLine[0:7] == "</head>":
                whereNow = 0
                continue
            
            if theLine[0:8] == "<article":
                whereNow = 2
                continue
            
            if theLine[0:10] == "</article>":
                whereNow = 0
                continue
            
            if whereNow == 1:
                if theLine[0:4] == "<!--":
                    htmlHead += theLine+"\n"
                    tmpBits   = theLine[4:-4].split("http://")
                    if len(tmpBits) == 2:
                        genStr = tmpBits[0].strip().replace("LaTeXML","<a href='http://"+tmpBits[1]+"'>LaTeXML</a>")
                if theLine[0:5] == "<link":
                    htmlHead += theLine+"\n"
                if theLine[0:7] == "<script":
                    htmlHead += theLine+"\n"
            
            if whereNow == 2:
                htmlBody += theLine+"\n"
    
    print("Done.")
    print("Writing PHP file ...")
    
    with open(path.join(htmlDir,phpFile),"w") as outFile:
        outFile.write("<?php\n")
        outFile.write("    $bMain  = true;\n")
        # outFile.write("    $sTitle = 'Manual';\n")
        outFile.write("    $sHead  = '"+htmlHead+"';\n")
        outFile.write("    $sHead .= '<link rel=\"stylesheet\" href=\"/SixTrack/css/latexml-fix.css\" type=\"text/css\">\n';\n")
        outFile.write("    require_once('../../includes/header.php');\n")
        outFile.write("?>\n")
        outFile.write("<article id='manual' class='ltx_document'>\n")
        outFile.write("<div id='generator'>"+genStr+"</div>\n")
        outFile.write(htmlBody)
        outFile.write("</article>\n")
        outFile.write("<aside>\n")
        outFile.write("    <?php include_once($incPath.'/aside_manual.php'); ?>\n")
        outFile.write("</aside>\n")
        outFile.write("<?php\n")
        outFile.write("    require_once($incPath.'/footer.php');\n")
        outFile.write("?>\n")
    
    print("Done.")
    
    return


if __name__ == "__main__":
    
    if len(sys.argv) == 2:
        htmlDir = sys.argv[1]
    else:
        print("Error: Input folder missing.")
        exit(1)
    
    htmlFiles = []
    for dirItem in listdir(htmlDir):
        if len(dirItem) < 5: continue
        if dirItem[-5:] == ".html":
            htmlFiles.append(dirItem[0:-5])
    # print(htmlFiles)
    #for htmlFile in htmlFiles:
    wrapHTML(htmlDir,"manual.html",htmlFiles)
    
    exit(0)
 
# End Main
