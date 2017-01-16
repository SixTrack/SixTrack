#This file performs the ast_mask filtering operations
#For the regex:
#sed -i -e '/^e  *'${OPT}'/d' ${CMAKE_BINARY_DIR}/ast_mask/${loop}.ast
FILE(READ ${FILE_NAME} INPUT_STRING)
STRING(REPLACE "e ${OPT}\n" "" OUTPUT_STRING ${INPUT_STRING} )
FILE(WRITE ${FILE_NAME} "${OUTPUT_STRING}")
