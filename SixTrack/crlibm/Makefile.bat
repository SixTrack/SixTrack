cd crlibm
fcc -c /K PENTIUM /O ericc.c
fcc -c /K PENTIUM /O logsix.c
fcc -c /K PENTIUM /O log_fast.c
fcc -c /K PENTIUM /O exp.c
fcc -c /K PENTIUM /O exp_fast.c
fcc -c /K PENTIUM /O atan.c
fcc -c /K PENTIUM /O atan_fast.c
fcc -c /K PENTIUM /O addition_scs.c
fcc -c /K PENTIUM /O multiplication_scs.c
fcc -c /K PENTIUM /O division_scs.c
fcc -c /K PENTIUM /O zero_scs.c
fcc -c /K PENTIUM /O double2scs.c
fcc -c /K PENTIUM /O scs2double.c
fcc -c /K PENTIUM /O sine.c
fcc -c /K PENTIUM /O cosine.c
fcc -c /K PENTIUM /O tan.c
fcc -c /K PENTIUM /O trigo_fast.c
fcc -c /K PENTIUM /O rem_pio2.c
fcc -c /K PENTIUM /O log10.c
fcc -c /K PENTIUM /O csh_fast.c
copy ericc.obj ..
copy logsix.obj ..
copy log_fast.obj ..
copy exp.obj ..
copy exp_fast.obj ..
copy atan.obj ..
copy atan_fast.obj ..
copy addition_scs.obj ..
copy multiplication_scs.obj ..
copy division_scs.obj ..
copy zero_scs.obj ..
copy double2scs.obj ..
copy scs2double.obj ..
copy sine.obj ..
copy cosine.obj ..
copy tan.obj ..
copy trigo_fast.obj ..
copy rem_pio2.obj ..
copy log10.obj ..
copy csh_fast.obj ..
gcc -c disable_xp.c
gcc -c enable_xp.c
gcc -m32 -std=c99 -W -Wall -pedantic -c round_near.c
gcc -m32 -std=c99 -W -Wall -pedantic -c dtostr.c
gcc -m32 -std=c99 -W -Wall -pedantic -c dtoa_c.c
copy disable_xp.o ..
copy enable_xp.o ..
copy round_near.o ..
copy dtoa_c.o ..
cd ..
