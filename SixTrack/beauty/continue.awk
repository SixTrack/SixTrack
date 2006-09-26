BEGIN {
  blanks="                                                                                         "
  nplines=0+0
  while (getline > 0) {
    l=$0
    if (l ~ /^+/) {
      nplines=nplines+1
      plines[nplines]=l
    }
    else {
      if (l ~ /^     &/) {
        if (comment == 1 || p == "") {
          print (FNR " ERROR!! Comment precedes CONTINUATION") 
        }
        p=substr(p blanks,1,72)"&"
      }
      print (p)
      for (i=1; i <= nplines; i++) {
        print (plines[i])
      }
      nplines=0
      p=l
      if (p ~ /^[^1-9 ]/) {
          comment=1
      }
      else {
        comment=0
      }
    }
  }
  print (p)
  for (i=1; i <= nplines; i++) {
    print (plines[i])
  }
}
