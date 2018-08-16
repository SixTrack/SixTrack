BEGIN {
  comment=0
  s=getline
  p=$0
  while (getline > 0) {
    l=$0
    if (l ~ /^+if crlibm/) {
      if (p ~ /^+if crlibm/) {
        print ("Line No: ",FNR)
      }
    }
    p=l
  }
}
