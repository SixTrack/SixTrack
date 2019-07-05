# Test Results

This test should produce the same results as the collimation_old-db_new-block test, but the results are not numerically identical.

The initial_state.dat file as well as the scatter files should be identical, otherwise the jobs are not configured correctly.

Expected differences in output files:

* The collgaps.dat file will have slightly different values of the angle due to conversion between radians and degrees.
* The final_state file and the dump files should only have very small differences in decimals. Since the coll DB is read slightly differently, some rounding differences are expected.

