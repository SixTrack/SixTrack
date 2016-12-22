BEGIN {
    extra_spaces = 0
    while (getline > 0) {
	#l = $0
	#print NF
	#print l
	if ( $0 ~ /^+/ ){ #Match astuce statement
	    #Remove comments
	    l = $0
	    sub(/!.*/, " ",l)

	    #Check that there are no extra spaces
	    split(l, a)
	    if ( length(a) > 2 ) {
		print NR ": '" l "'"
	        extra_spaces = extra_spaces + 1
	    }
	}
    }
    if (extra_spaces != 0){
	print ("Error with +if's - there should be no spaces in their arguments!")
	exit 1
    }
}
