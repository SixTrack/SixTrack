BEGIN {
  while (getline > 0) {
    l=$0
    if (l !~ /\/\//) {
      print (l)
    }
    else {
      len=length(l)
      newl=""
      n=split(l,strings,"\"")     
      for (i=1; i <= n; i++) {
        if (i%2 == 0) {
          newl=newl "\"" (strings[i]) "\""
        }
        else {
          newstring=""
          stat=match(strings[i],"//")
          if (stat == 0) {
            newstring=strings[i]
          }
          else {
            ll=length(strings[i])
            newstring=substr(strings[i],1,RSTART-1) 
            newstring=newstring " /* "
            newstring=newstring substr(strings[i],RSTART+2,ll) " */ "
          }
          newl=newl newstring
        }
      }
      print (newl)
    }
  }
}
