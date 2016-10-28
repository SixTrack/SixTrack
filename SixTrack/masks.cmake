#This file performs the ast_mask filtering operations
#For the regex:
#sed -e 's/\(^df .*\), *'$opt' *,\(.*\)/\1,\2/'
#	-e 's/\(^df .*\), *'$opt' *$/\1/'
#	-e 's/\(^df  *'$opt' *,\)\(.*\)/df \2/'
#	-e 's/\(^df  *'$opt' *$\)/df/'
#First entry removes ,OPT, (middle of a line) 's/\(^df .*\), *${OPT} *,\(.*\)/,/g'
#Second entry removes ,OPT (last entry on a line)
#Third entry removes OPT, (first entry on a line)
#Last entry removes OPT (alone on a line) 's/\(^df  *'${OPT}' *$\)/df/'  \
FILE(READ ${FILE_NAME} INPUT_STRING)
STRING(REPLACE ",${OPT}," "," TEMP_A ${INPUT_STRING} )
STRING(REPLACE ",${OPT}\n" "\n" TEMP_B ${TEMP_A} )
STRING(REPLACE " ${OPT}," " " OUTPUT_STRING ${TEMP_B} )
FILE(WRITE ${FILE_NAME} "${OUTPUT_STRING}")
