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
lf95.exe -c -o1 -tp -fix -noconcc asin_rn.f
lf95.exe -c -o1 -tp -fix -noconcc acos_rn.f
lf95.exe -c -o1 -tp -fix -noconcc atan2_rn.f
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
copy asin_rn.obj ..
copy acos_rn.obj ..
copy atan2_rn.obj ..
cd ..
lf95.exe -c -o1 -tp -fix -nconcc sixve.f
lf95.exe -c -o1 -tp -fix -nconcc sixvefox.f
lf95.exe -c -o1 -tp -fix -nconcc track.f
lf95.exe -c -o1 -tp -fix -nconcc dabnews.f
lf95.exe -c -o1 -tp -fix -nconcc lielib.f
lf95.exe -c -o1 -tp -fix -nconcc myboinc.f
lf95.exe -out %%SIXTRACK%%.exe %%BOINC%% sixve.obj sixvefox.obj track.obj dabnews.obj lielib.obj ericc.obj atan2_rn.obj atan.obj atan_fast.obj exp.obj logsix.obj exp_fast.obj log_fast.obj addition_scs.obj multiplication_scs.obj zero_scs.obj double2scs.obj scs2double.obj division_scs.obj sine.obj cosine.obj tan.obj trigo_fast.obj rem_pio2.obj log10.obj asin_rn.obj acos_rn.obj csh_fast.obj
