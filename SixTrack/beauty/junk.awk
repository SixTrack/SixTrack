BEGIN {
  plus=0
  comment=0
  blanks="                                                                                         "
  s=getline
  prep=$0
  s=getline
  p=$0
  while (getline > 0) {
    l=$0
    if (l == "" || l !~ /^[0-9 +]/) {
      print (prep)
      if (l !~ /^+/) {
        comment=1
      }
      else {
        plus=1
      }
    }
    else {
      if (l ~ /^     &/) {
        if (plus == 1) {
          prep=substr(prep blanks,1,72)"&"
        }
        else { 
          p=substr(p blanks,1,72)"&"
        }
      }
      else {
      }
      comment=0
      plus=0
      print(prep)
    }
    prep=p
    p=l
  }
  print (prep)
  print (p)
}
