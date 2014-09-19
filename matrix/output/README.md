Each subdirectory contains root files for the following tests

Note: in all tests, bXings with `h_b * h_y = 0` are omitted, as well as those 
considered deviant ("kicked") from the rellum analysis

** EVERYTHING HERE IS FOR RUN 13 **

- `all_bXings`: no additional bXings omitted from bunch fitting

- `spinbit_equalised`: randomly removed the minimum number of bXings such that the 
  number of bXings per spin state is equal

- `24_bXings`: only consider bXings 72-103 (4 repetitions of spin pattern; 72 = 0 (mod 8) )
  - see log entry 14.09.12

- `omit_8_bXings_after_aborts_69_70`: omitted the first 8 bXings after the abort gaps as well
  as 8 bXings after bXings 69 & 70 (which were empty for the first ~230 runs)
  - see log entry 14.09.16

- `omit_8_bXings_before_aborts_69_70`: omitted the last 8 bXings before the abort gaps as well 
  as 8 bXings before bXings 69 & 70
  - see log entry 14.09.16

- `omit_2_bXings_after_aborts_69_70`: omitted first 2 bXings after abort gaps as well as 
  first 2 bXings after bXings 69 & 70
  - see log entry 14.09.17
