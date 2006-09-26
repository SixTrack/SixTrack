BEGIN {
  count=0+0
  mcount=0+0
  nested=0+0
  indent=""
  p=""
  pif=""
  while (getline > 0 && count >= 0) {
    l=$0
    if (l ~ /^+if/) {
      count=count+1
      if (count == 1) {
        pif=l
        pfnr=FNR
      }
      if (count == 2 && nested == 0) {
        print ("") > "blah"
        print (pfnr) > "blah"
        print (pif) > "blah"
        nested=1
      }
      if (count >= 2) {
        indent = indent"     "
        print (indent FNR) > "blah"
        print (indent l) > "blah"
      }
      if (count > mcount) {
        mcount=count
      }
    }
    if (l ~ /^+ei/) {
      count=count-1
      if (count >= 1) {
        print (indent FNR) > "blah"
        print (indent l) > "blah"
        len=length(indent)-5
        indent=substr(indent,1,len)
      }
      if (count < 0) {
        print ("ERROR : Count negative !!!")
        print ("Line No: " FNR l)
      }
      if (count == 0) {
        pif=""
        if (nested == 1) {
          nested=0
          print (indent FNR) > "blah"
          print (indent l) > "blah"
        }
      }
    }
  p=l
  }
  if (count != 0) {
    print ("Error!!! Count=" count)
  }
  print ("Maximum nesting level=" mcount)
}
