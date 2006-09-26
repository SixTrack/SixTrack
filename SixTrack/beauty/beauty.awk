BEGIN {
  while (getline > 0) {
    l=$0
    if (l ~ /^[^1-9 ]/) {
      print (l)
    }
    else {
      len=length(l)
      newl=""
      n=split(l,strings,"'")     
      for (i=1; i <= n; i++) {
        if (i%2 == 0) {
          newl=newl "'" (strings[i]) "'"
        }
        else {
          newl=newl tolower(strings[i])
        }
      }
      print (newl)
    }
  }
}
