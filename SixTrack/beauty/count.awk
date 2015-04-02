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
        print ("") > "count.log"
        print (pfnr) > "count.log"
        print (pif) > "count.log"
        nested=1
      }
      if (count >= 2) {
        indent = indent"     "
        print (indent FNR) > "count.log"
        print (indent l) > "count.log"
      }
      if (count > mcount) {
        mcount=count
      }
    }
    if (l ~ /^+ei/) {
      count=count-1
      if (count >= 1) {
        print (indent FNR) > "count.log"
        print (indent l) > "count.log"
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
          print (indent FNR) > "count.log"
          print (indent l) > "count.log"
        }
      }
    }
  p=l
  }
  if (count != 0) {
    print ("Error!!! Count=" count)
    exit 1
  }
  print ("Maximum nesting level=" mcount)
}
