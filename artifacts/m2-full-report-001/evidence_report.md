# M2 hypothesis evidence

Endpoint parameter hash: da822ce1b48d242f1585fbf4d1811eec
Archived Phase 27 code MD5: d168daef9590c4a0874ba76ca163f1ee
Source scope: HISTORICAL_CODE

## Coverage
 season_id policy h1_matched_beta_only                h2_matched_S_beta
   2013/14    raw            COMPLETED                        COMPLETED
   2013/14  blend            COMPLETED                        COMPLETED
   2014/15    raw            COMPLETED                        COMPLETED
   2014/15  blend            COMPLETED                        COMPLETED
   2015/16    raw            COMPLETED                        COMPLETED
   2015/16  blend            COMPLETED                        COMPLETED
   2016/17    raw            COMPLETED                        COMPLETED
   2016/17  blend            COMPLETED                        COMPLETED
   2017/18    raw            COMPLETED                        COMPLETED
   2017/18  blend            COMPLETED                        COMPLETED
   2018/19    raw            COMPLETED                        COMPLETED
   2018/19  blend            COMPLETED                        COMPLETED
   2019/20    raw            COMPLETED                        COMPLETED
   2019/20  blend            COMPLETED                        COMPLETED
   2020/21    raw            COMPLETED COMPLETED_NUMERICALLY_INELIGIBLE
   2020/21  blend            COMPLETED COMPLETED_NUMERICALLY_INELIGIBLE
   2021/22    raw            COMPLETED                        COMPLETED
   2021/22  blend            COMPLETED                        COMPLETED
   2022/23    raw            COMPLETED COMPLETED_NUMERICALLY_INELIGIBLE
   2022/23  blend            COMPLETED COMPLETED_NUMERICALLY_INELIGIBLE
 h3_local_I                 h3_finite_I_path
  COMPLETED                        COMPLETED
  COMPLETED                        COMPLETED
  COMPLETED                        COMPLETED
  COMPLETED                        COMPLETED
  COMPLETED                        COMPLETED
  COMPLETED                        COMPLETED
  COMPLETED                        COMPLETED
  COMPLETED                        COMPLETED
  COMPLETED                        COMPLETED
  COMPLETED                        COMPLETED
  COMPLETED                        COMPLETED
  COMPLETED                        COMPLETED
  COMPLETED                        COMPLETED
  COMPLETED                        COMPLETED
  COMPLETED COMPLETED_NUMERICALLY_INELIGIBLE
  COMPLETED COMPLETED_NUMERICALLY_INELIGIBLE
  COMPLETED COMPLETED_NUMERICALLY_INELIGIBLE
  COMPLETED COMPLETED_NUMERICALLY_INELIGIBLE
  COMPLETED COMPLETED_NUMERICALLY_INELIGIBLE
  COMPLETED COMPLETED_NUMERICALLY_INELIGIBLE

## H1
Saved refined intervals: 6 dominance exchanges, 173 displacements.
Original-grid timing-event bounds by policy and path:
 policy              path events resolved_displacement possible_displacement
    raw MATCHED_BETA_ONLY     80                    72                    79
  blend MATCHED_BETA_ONLY     77                    69                    75
 unresolved mixed
          5     4
          5     3
Targeted signed competitor-height reversals:
 policy season              path                origin   actual_window
    raw   2017 MATCHED_BETA_ONLY MATCHED_SEASON_ORIGIN EPIDEMIC_SEASON
    raw   2018 MATCHED_BETA_ONLY MATCHED_SEASON_ORIGIN EPIDEMIC_SEASON
    raw   2022 MATCHED_BETA_ONLY MATCHED_SEASON_ORIGIN EPIDEMIC_SEASON
  blend   2017 MATCHED_BETA_ONLY MATCHED_SEASON_ORIGIN EPIDEMIC_SEASON
  blend   2018 MATCHED_BETA_ONLY MATCHED_SEASON_ORIGIN EPIDEMIC_SEASON
  blend   2022 MATCHED_BETA_ONLY MATCHED_SEASON_ORIGIN EPIDEMIC_SEASON
      v_left     v_right left_week right_week early_minus_late_left
 -0.16059965 -0.16024356        24         11             -4.485494
 -0.16914598 -0.16878988        23         12             -4.297457
  0.06694620  0.06837058        51          4            -38.573675
 -0.16309233 -0.16273623        24         11             -6.666722
 -0.16166794 -0.16131185        22         12             -4.479822
  0.06837058  0.06979497        51          4            -29.218166
 early_minus_late_right left_branch right_branch
              2.5047931   tracked-3    tracked-2
              4.7285323   tracked-2    tracked-1
              3.0208547   tracked-3    tracked-1
              0.3482759   tracked-2    tracked-1
              4.2390798   tracked-3    tracked-2
             13.3642483   tracked-3    tracked-1
The bounds keep mixed and unresolved parent intervals. Blending is a numerical robustness arm.
 season_id policy shift_events displacement_events exchange_events mixed_events
   2013/14    raw           14                  14               0            0
   2013/14  blend           13                  13               0            0
   2014/15    raw           12                  12               0            0
   2014/15  blend           12                  12               0            0
   2015/16    raw           13                  13               0            0
   2015/16  blend           13                  13               0            0
   2016/17    raw           13                  13               0            0
   2016/17  blend           12                  12               0            0
   2017/18    raw            3                   2               0            1
   2017/18  blend            3                   2               1            0
   2018/19    raw            4                   3               0            1
   2018/19  blend            4                   3               0            1
   2019/20    raw            6                   5               0            1
   2019/20  blend            5                   4               0            1
   2020/21    raw           12                  10               0            1
   2020/21  blend           12                  10               0            1
   2021/22    raw            1                   0               0            0
   2021/22  blend            1                   0               0            0
   2022/23    raw            2                   0               1            0
   2022/23  blend            2                   0               1            0
 unresolved_events no_shift displacement_movement exchange_movement
                 0    FALSE                    15                 0
                 0    FALSE                    15                 0
                 0    FALSE                    16                 0
                 0    FALSE                    17                 0
                 0    FALSE                    14                 0
                 0    FALSE                    14                 0
                 0    FALSE                    14                 0
                 0    FALSE                    14                 0
                 0    FALSE                     2                 0
                 0    FALSE                     2                13
                 0    FALSE                     3                 0
                 0    FALSE                     3                 0
                 1    FALSE                     5                 0
                 1    FALSE                     4                 0
                 2    FALSE                    13                 0
                 2    FALSE                    13                 0
                 1    FALSE                     0                 0
                 1    FALSE                     0                 0
                 1    FALSE                     0                47
                 1    FALSE                     0                47
 mixed_parent_movement unresolved_parent_movement total_absolute_movement
                     0                          0                      15
                     0                          0                      15
                     0                          0                      16
                     0                          0                      17
                     0                          0                      14
                     0                          0                      14
                     0                          0                      14
                     0                          0                      14
                    14                          0                      16
                     0                          0                      15
                    12                          0                      15
                    11                          0                      14
                     7                          7                      12
                     7                          7                      11
                     4                          5                      18
                     3                          4                      17
                     0                          1                       1
                     0                          1                       1
                     0                          1                      48
                     0                          1                      48
 signed_endpoint_movement identified_domain unresolved_domain
                      -15         0.3646431       0.000000000
                      -15         0.3646431       0.000000000
                      -16         0.3646431       0.000000000
                      -17         0.3646431       0.000000000
                      -14         0.3646431       0.000000000
                      -14         0.3646431       0.000000000
                      -14         0.3646431       0.000000000
                      -14         0.3646431       0.000000000
                      -16         0.3646431       0.000000000
                      -15         0.3646431       0.000000000
                      -15         0.3646431       0.000000000
                      -14         0.3646431       0.000000000
                      -12         0.3190627       0.045580389
                      -11         0.3190627       0.045580389
                      -18         0.3162140       0.048429164
                      -17         0.3162140       0.048429164
                        1         0.3589456       0.005697549
                        1         0.3589456       0.005697549
                      -48         0.3617943       0.002848774
                      -48         0.3617943       0.002848774
 displacement_event_fraction
                   1.0000000
                   1.0000000
                   1.0000000
                   1.0000000
                   1.0000000
                   1.0000000
                   1.0000000
                   1.0000000
                   1.0000000
                   0.6666667
                   1.0000000
                   1.0000000
                   1.0000000
                   1.0000000
                   1.0000000
                   1.0000000
                          NA
                          NA
                   0.0000000
                   0.0000000
 policy season              path          classification absolute_movement
  blend   2020 MATCHED_BETA_ONLY      BOUNDARY_CANDIDATE                 1
    raw   2020 MATCHED_BETA_ONLY      BOUNDARY_CANDIDATE                 1
  blend   2021 MATCHED_BETA_ONLY      BOUNDARY_CANDIDATE                 1
    raw   2021 MATCHED_BETA_ONLY      BOUNDARY_CANDIDATE                 1
  blend   2022 MATCHED_BETA_ONLY      BOUNDARY_CANDIDATE                 1
    raw   2022 MATCHED_BETA_ONLY      BOUNDARY_CANDIDATE                 1
  blend   2013 MATCHED_BETA_ONLY CONTINUOUS_DISPLACEMENT                15
    raw   2013 MATCHED_BETA_ONLY CONTINUOUS_DISPLACEMENT                15
  blend   2014 MATCHED_BETA_ONLY CONTINUOUS_DISPLACEMENT                17
    raw   2014 MATCHED_BETA_ONLY CONTINUOUS_DISPLACEMENT                16
  blend   2015 MATCHED_BETA_ONLY CONTINUOUS_DISPLACEMENT                14
    raw   2015 MATCHED_BETA_ONLY CONTINUOUS_DISPLACEMENT                14
  blend   2016 MATCHED_BETA_ONLY CONTINUOUS_DISPLACEMENT                14
    raw   2016 MATCHED_BETA_ONLY CONTINUOUS_DISPLACEMENT                14
  blend   2017 MATCHED_BETA_ONLY CONTINUOUS_DISPLACEMENT                 2
    raw   2017 MATCHED_BETA_ONLY CONTINUOUS_DISPLACEMENT                 3
  blend   2018 MATCHED_BETA_ONLY CONTINUOUS_DISPLACEMENT                 4
    raw   2018 MATCHED_BETA_ONLY CONTINUOUS_DISPLACEMENT                 4
  blend   2019 MATCHED_BETA_ONLY CONTINUOUS_DISPLACEMENT                 8
    raw   2019 MATCHED_BETA_ONLY CONTINUOUS_DISPLACEMENT                 7
  blend   2020 MATCHED_BETA_ONLY CONTINUOUS_DISPLACEMENT                15
    raw   2020 MATCHED_BETA_ONLY CONTINUOUS_DISPLACEMENT                15
  blend   2017 MATCHED_BETA_ONLY      DOMINANCE_EXCHANGE                13
    raw   2017 MATCHED_BETA_ONLY      DOMINANCE_EXCHANGE                13
  blend   2018 MATCHED_BETA_ONLY      DOMINANCE_EXCHANGE                10
    raw   2018 MATCHED_BETA_ONLY      DOMINANCE_EXCHANGE                11
  blend   2022 MATCHED_BETA_ONLY      DOMINANCE_EXCHANGE                47
    raw   2022 MATCHED_BETA_ONLY      DOMINANCE_EXCHANGE                47
  blend   2013 MATCHED_BETA_ONLY      SAMPLED_PEAK_BIRTH                 0
    raw   2013 MATCHED_BETA_ONLY      SAMPLED_PEAK_BIRTH                 0
  blend   2014 MATCHED_BETA_ONLY      SAMPLED_PEAK_BIRTH                 0
  blend   2016 MATCHED_BETA_ONLY      SAMPLED_PEAK_BIRTH                 0
    raw   2016 MATCHED_BETA_ONLY      SAMPLED_PEAK_BIRTH                 0
  blend   2017 MATCHED_BETA_ONLY      SAMPLED_PEAK_BIRTH                 0
    raw   2017 MATCHED_BETA_ONLY      SAMPLED_PEAK_BIRTH                 0
  blend   2018 MATCHED_BETA_ONLY      SAMPLED_PEAK_BIRTH                 0
    raw   2018 MATCHED_BETA_ONLY      SAMPLED_PEAK_BIRTH                 0
  blend   2019 MATCHED_BETA_ONLY      SAMPLED_PEAK_BIRTH                 0
    raw   2019 MATCHED_BETA_ONLY      SAMPLED_PEAK_BIRTH                 0
  blend   2020 MATCHED_BETA_ONLY      SAMPLED_PEAK_BIRTH                 1
    raw   2020 MATCHED_BETA_ONLY      SAMPLED_PEAK_BIRTH                 2
    raw   2021 MATCHED_BETA_ONLY      SAMPLED_PEAK_BIRTH                 0
  blend   2013 MATCHED_BETA_ONLY       SAMPLED_PEAK_LOSS                 0
  blend   2016 MATCHED_BETA_ONLY       SAMPLED_PEAK_LOSS                 0
    raw   2016 MATCHED_BETA_ONLY       SAMPLED_PEAK_LOSS                 0
  blend   2017 MATCHED_BETA_ONLY       SAMPLED_PEAK_LOSS                 0
    raw   2017 MATCHED_BETA_ONLY       SAMPLED_PEAK_LOSS                 0
  blend   2018 MATCHED_BETA_ONLY       SAMPLED_PEAK_LOSS                 0
    raw   2018 MATCHED_BETA_ONLY       SAMPLED_PEAK_LOSS                 0
  blend   2019 MATCHED_BETA_ONLY       SAMPLED_PEAK_LOSS                 0
    raw   2019 MATCHED_BETA_ONLY       SAMPLED_PEAK_LOSS                 0
  blend   2013 MATCHED_BETA_ONLY      STABLE_SAME_BRANCH                 0
    raw   2013 MATCHED_BETA_ONLY      STABLE_SAME_BRANCH                 0
  blend   2014 MATCHED_BETA_ONLY      STABLE_SAME_BRANCH                 0
    raw   2014 MATCHED_BETA_ONLY      STABLE_SAME_BRANCH                 0
  blend   2015 MATCHED_BETA_ONLY      STABLE_SAME_BRANCH                 0
    raw   2015 MATCHED_BETA_ONLY      STABLE_SAME_BRANCH                 0
  blend   2016 MATCHED_BETA_ONLY      STABLE_SAME_BRANCH                 0
    raw   2016 MATCHED_BETA_ONLY      STABLE_SAME_BRANCH                 0
  blend   2017 MATCHED_BETA_ONLY      STABLE_SAME_BRANCH                 0
    raw   2017 MATCHED_BETA_ONLY      STABLE_SAME_BRANCH                 0
  blend   2018 MATCHED_BETA_ONLY      STABLE_SAME_BRANCH                 0
    raw   2018 MATCHED_BETA_ONLY      STABLE_SAME_BRANCH                 0
  blend   2019 MATCHED_BETA_ONLY      STABLE_SAME_BRANCH                 0
    raw   2019 MATCHED_BETA_ONLY      STABLE_SAME_BRANCH                 0
  blend   2020 MATCHED_BETA_ONLY      STABLE_SAME_BRANCH                 0
    raw   2020 MATCHED_BETA_ONLY      STABLE_SAME_BRANCH                 0
  blend   2021 MATCHED_BETA_ONLY      STABLE_SAME_BRANCH                 0
    raw   2021 MATCHED_BETA_ONLY      STABLE_SAME_BRANCH                 0
  blend   2022 MATCHED_BETA_ONLY      STABLE_SAME_BRANCH                 0
    raw   2022 MATCHED_BETA_ONLY      STABLE_SAME_BRANCH                 0
  blend   2019 MATCHED_BETA_ONLY              UNRESOLVED                 3
    raw   2019 MATCHED_BETA_ONLY              UNRESOLVED                 5
Absolute movement may exceed signed net endpoint movement. Equal-season summaries retain no-shift seasons.

## H2
 reference policy                design attempted eligible beta_larger
        M2    raw           2013_ORIGIN        55       49          49
        M2    raw MATCHED_SEASON_ORIGIN        50       39          35
        M2  blend           2013_ORIGIN        55       49          49
        M2  blend MATCHED_SEASON_ORIGIN        50       39          35
     JOINT    raw           2013_ORIGIN        55       42          42
     JOINT    raw MATCHED_SEASON_ORIGIN        50       40          35
 negative_S negative_beta
         19            44
         35            35
         19            44
         35            35
         42            42
         35            35
 season_id policy tested_steps eligible_steps  selected_h selected_E_S
   2013/14    raw            5              5 0.011395097   -25.431004
   2013/14  blend            5              5 0.011395097   -27.691637
   2014/15    raw            5              5 0.011395097   -23.614165
   2014/15  blend            5              5 0.011395097   -24.459497
   2015/16    raw            5              5 0.011395097   -28.241045
   2015/16  blend            5              5 0.011395097   -32.110992
   2016/17    raw            5              5 0.011395097   -24.771859
   2016/17  blend            5              5 0.011395097   -24.254684
   2017/18    raw            5              5 0.011395097    -6.673089
   2017/18  blend            5              5 0.011395097    -6.718536
   2018/19    raw            5              5 0.011395097    -7.584221
   2018/19  blend            5              5 0.011395097    -7.882479
   2019/20    raw            5              5 0.011395097    -9.784061
   2019/20  blend            5              5 0.011395097   -11.081076
   2020/21    raw            5              0          NA           NA
   2020/21  blend            5              0          NA           NA
   2021/22    raw            5              4 0.005697549    11.101892
   2021/22  blend            5              4 0.005697549     9.623076
   2022/23    raw            5              0          NA           NA
   2022/23  blend            5              0          NA           NA
 selected_E_beta magnitude_order    direction           selection_rule
      -31.243127     BETA_LARGER BOTH_ADVANCE MEDIAN_ELIGIBLE_LOG_STEP
      -33.482601     BETA_LARGER BOTH_ADVANCE MEDIAN_ELIGIBLE_LOG_STEP
      -29.685055     BETA_LARGER BOTH_ADVANCE MEDIAN_ELIGIBLE_LOG_STEP
      -30.215301     BETA_LARGER BOTH_ADVANCE MEDIAN_ELIGIBLE_LOG_STEP
      -34.127620     BETA_LARGER BOTH_ADVANCE MEDIAN_ELIGIBLE_LOG_STEP
      -38.457075     BETA_LARGER BOTH_ADVANCE MEDIAN_ELIGIBLE_LOG_STEP
      -29.121108     BETA_LARGER BOTH_ADVANCE MEDIAN_ELIGIBLE_LOG_STEP
      -28.576873     BETA_LARGER BOTH_ADVANCE MEDIAN_ELIGIBLE_LOG_STEP
       -8.497160     BETA_LARGER BOTH_ADVANCE MEDIAN_ELIGIBLE_LOG_STEP
       -8.541691     BETA_LARGER BOTH_ADVANCE MEDIAN_ELIGIBLE_LOG_STEP
       -9.510873     BETA_LARGER BOTH_ADVANCE MEDIAN_ELIGIBLE_LOG_STEP
       -9.846208     BETA_LARGER BOTH_ADVANCE MEDIAN_ELIGIBLE_LOG_STEP
      -11.985708     BETA_LARGER BOTH_ADVANCE MEDIAN_ELIGIBLE_LOG_STEP
      -13.432761     BETA_LARGER BOTH_ADVANCE MEDIAN_ELIGIBLE_LOG_STEP
              NA      INELIGIBLE   INELIGIBLE MEDIAN_ELIGIBLE_LOG_STEP
              NA      INELIGIBLE   INELIGIBLE MEDIAN_ELIGIBLE_LOG_STEP
        7.020422        S_LARGER   BOTH_DELAY MEDIAN_ELIGIBLE_LOG_STEP
        5.683791        S_LARGER   BOTH_DELAY MEDIAN_ELIGIBLE_LOG_STEP
              NA      INELIGIBLE   INELIGIBLE MEDIAN_ELIGIBLE_LOG_STEP
              NA      INELIGIBLE   INELIGIBLE MEDIAN_ELIGIBLE_LOG_STEP
 step_range_S step_range_beta
   1.27313627      2.02850426
   2.24686117      4.65920405
   0.55457258      0.84081610
   0.60613245      1.04491414
   1.25204191      1.59257246
   1.46444348      2.22428280
  11.45312650     21.77171264
   7.63699990     14.18538079
   0.10101188      0.15167310
   0.08109180      0.14422148
   0.06694466      0.14372448
   0.16945466      0.26424345
   0.22229265      0.24766424
   0.49123070      0.67892778
           NA              NA
           NA              NA
   0.09836972      0.09925816
   0.07106313      0.07035684
           NA              NA
           NA              NA
Eligible matched-origin magnitude-order exceptions:
 policy                design year           h       E_S   E_beta
    raw MATCHED_SEASON_ORIGIN 2021 0.022790195 11.195622 7.115005
    raw MATCHED_SEASON_ORIGIN 2021 0.011395097 11.120492 7.039172
    raw MATCHED_SEASON_ORIGIN 2021 0.005697549 11.101892 7.020422
    raw MATCHED_SEASON_ORIGIN 2021 0.002848774 11.097253 7.015747
  blend MATCHED_SEASON_ORIGIN 2021 0.022790195  9.690783 5.750828
  blend MATCHED_SEASON_ORIGIN 2021 0.011395097  9.636528 5.697100
  blend MATCHED_SEASON_ORIGIN 2021 0.005697549  9.623076 5.683791
  blend MATCHED_SEASON_ORIGIN 2021 0.002848774  9.619720 5.680471
Magnitude and sign are separate conclusions. Rows are year-step diagnostics, not independent samples.
Matched-origin 2021/22 contains the saved magnitude-order exceptions; 2013-origin 2017 slopes can be positive.

## H3
State-memory scenarios supplied: 1514
 season_id policy tested_steps eligible_steps selected_h selected_E_I
   2013/14    raw            5              5  0.0113951    -1.165382
   2013/14  blend            5              5  0.0113951    -1.079800
   2014/15    raw            5              5  0.0113951    -1.300696
   2014/15  blend            5              5  0.0113951    -1.230674
   2015/16    raw            5              5  0.0113951    -2.483109
   2015/16  blend            5              5  0.0113951    -2.877719
   2016/17    raw            5              5  0.0113951    -5.389350
   2016/17  blend            5              5  0.0113951    -5.169527
   2017/18    raw            5              5  0.0113951    -2.559195
   2017/18  blend            5              5  0.0113951    -2.541851
   2018/19    raw            5              5  0.0113951    -2.472416
   2018/19  blend            5              5  0.0113951    -2.492809
   2019/20    raw            5              5  0.0113951    -2.494985
   2019/20  blend            5              5  0.0113951    -2.699165
   2020/21    raw            5              0         NA           NA
   2020/21  blend            5              0         NA           NA
   2021/22    raw            5              5  0.0113951    -4.682854
   2021/22  blend            5              5  0.0113951    -4.500835
   2022/23    raw            5              0         NA           NA
   2022/23  blend            5              0         NA           NA
           selection_rule step_range_I                    reason
 MEDIAN_ELIGIBLE_LOG_STEP 4.315718e-04                          
 MEDIAN_ELIGIBLE_LOG_STEP 2.501232e-04                          
 MEDIAN_ELIGIBLE_LOG_STEP 4.580634e-04                          
 MEDIAN_ELIGIBLE_LOG_STEP 8.671261e-05                          
 MEDIAN_ELIGIBLE_LOG_STEP 2.652178e-03                          
 MEDIAN_ELIGIBLE_LOG_STEP 2.198488e-01                          
 MEDIAN_ELIGIBLE_LOG_STEP 1.259969e-02                          
 MEDIAN_ELIGIBLE_LOG_STEP 1.189330e-02                          
 MEDIAN_ELIGIBLE_LOG_STEP 3.144687e-03                          
 MEDIAN_ELIGIBLE_LOG_STEP 3.134641e-03                          
 MEDIAN_ELIGIBLE_LOG_STEP 3.521289e-03                          
 MEDIAN_ELIGIBLE_LOG_STEP 3.477628e-03                          
 MEDIAN_ELIGIBLE_LOG_STEP 4.052438e-03                          
 MEDIAN_ELIGIBLE_LOG_STEP 3.370906e-03                          
 MEDIAN_ELIGIBLE_LOG_STEP           NA PEAK_OR_BRANCH_INELIGIBLE
 MEDIAN_ELIGIBLE_LOG_STEP           NA PEAK_OR_BRANCH_INELIGIBLE
 MEDIAN_ELIGIBLE_LOG_STEP 3.397602e-02                          
 MEDIAN_ELIGIBLE_LOG_STEP 3.024816e-02                          
 MEDIAN_ELIGIBLE_LOG_STEP           NA PEAK_OR_BRANCH_INELIGIBLE
 MEDIAN_ELIGIBLE_LOG_STEP           NA PEAK_OR_BRANCH_INELIGIBLE
 policy year      low_z    high_z I_low_ratio I_high_ratio baseline_time
    raw 2013 -0.1823216 0.1823216   0.8333333          1.2     21.373228
    raw 2014 -0.1823216 0.1823216   0.8333333          1.2     23.144973
    raw 2015 -0.1823216 0.1823216   0.8333333          1.2     21.367486
    raw 2016 -0.1823216 0.1823216   0.8333333          1.2     14.065668
    raw 2017 -0.1823216 0.1823216   0.8333333          1.2     10.302921
    raw 2018 -0.1823216 0.1823216   0.8333333          1.2     10.066259
    raw 2019 -0.1823216 0.1823216   0.8333333          1.2     10.337169
    raw 2020 -0.1823216 0.1823216   0.8333333          1.2            NA
    raw 2021 -0.1823216 0.1823216   0.8333333          1.2      1.893794
    raw 2022 -0.1823216 0.1823216   0.8333333          1.2            NA
  blend 2013 -0.1823216 0.1823216   0.8333333          1.2     21.169605
  blend 2014 -0.1823216 0.1823216   0.8333333          1.2     22.704919
  blend 2015 -0.1823216 0.1823216   0.8333333          1.2     20.510971
  blend 2016 -0.1823216 0.1823216   0.8333333          1.2     14.021451
  blend 2017 -0.1823216 0.1823216   0.8333333          1.2     10.275595
  blend 2018 -0.1823216 0.1823216   0.8333333          1.2     10.239890
  blend 2019 -0.1823216 0.1823216   0.8333333          1.2     10.623574
  blend 2020 -0.1823216 0.1823216   0.8333333          1.2            NA
  blend 2021 -0.1823216 0.1823216   0.8333333          1.2      1.989126
  blend 2022 -0.1823216 0.1823216   0.8333333          1.2            NA
 low_to_baseline baseline_to_high timing_change discrete_change centroid_change
      -0.2122954       -0.2050562    -0.4173516              -1     -0.42758730
      -0.2469665       -0.2300234    -0.4769899               0     -0.64603775
      -0.4243344       -0.4544357    -0.8787701              -1     -0.97745076
      -1.3020918       -0.8417211    -2.1438129              -2     -1.18013945
      -0.5020820       -0.4618826    -0.9639646              -1     -0.52933973
      -0.4641734       -0.4611751    -0.9253485              -1      0.02644887
      -0.4808714       -0.4466752    -0.9275466              -1     -0.34106788
              NA               NA            NA              NA     -0.65558880
              NA               NA            NA              NA      0.27585615
              NA               NA            NA              NA      0.18164910
      -0.2034833       -0.1917388    -0.3952221               0     -0.36781995
      -0.2280382       -0.2221558    -0.4501940              -1     -0.59909196
      -0.4669379       -0.5998146    -1.0667525              -1     -0.97257508
      -1.1934111       -0.8069857    -2.0003968              -2     -1.17396120
      -0.4977401       -0.4602643    -0.9580044              -1     -0.53886919
      -0.4868359       -0.4525651    -0.9394010              -1      0.02416873
      -0.4684525       -0.4655136    -0.9339661              -1     -0.35673117
              NA               NA            NA              NA     -0.62325574
              NA               NA            NA              NA      0.21602352
              NA               NA            NA              NA      0.20058468
 height_change burden_change                peak_status at_least_one_week
    -782.42048     -729.9966                   ELIGIBLE             FALSE
    -760.48744    -1197.4135                   ELIGIBLE             FALSE
    -537.16880    -1431.5367                   ELIGIBLE             FALSE
     332.21494    -1883.0920                   ELIGIBLE              TRUE
     370.24402    -1072.8702                   ELIGIBLE             FALSE
     237.82923      655.4166                   ELIGIBLE             FALSE
     -47.02143     -471.8655                   ELIGIBLE             FALSE
    -959.11545     4273.1947 INELIGIBLE_PEAK_OR_SUPPORT                NA
    3057.22699      761.8673 INELIGIBLE_PEAK_OR_SUPPORT                NA
   -1635.12094     2517.6668 INELIGIBLE_PEAK_OR_SUPPORT                NA
    -768.21341     -602.4418                   ELIGIBLE             FALSE
    -850.45566    -1128.0786                   ELIGIBLE             FALSE
    -570.32357    -1505.4195                   ELIGIBLE              TRUE
     340.13425    -1865.0292                   ELIGIBLE              TRUE
     369.54782    -1093.5596                   ELIGIBLE             FALSE
     240.62212      648.7385                   ELIGIBLE             FALSE
     -55.58701     -524.4600                   ELIGIBLE             FALSE
    -878.69435     4705.6671 INELIGIBLE_PEAK_OR_SUPPORT                NA
    3056.34287      636.9746 INELIGIBLE_PEAK_OR_SUPPORT                NA
   -1649.56682     2611.7348 INELIGIBLE_PEAK_OR_SUPPORT                NA
 at_least_two_weeks at_least_four_weeks
              FALSE               FALSE
              FALSE               FALSE
              FALSE               FALSE
               TRUE               FALSE
              FALSE               FALSE
              FALSE               FALSE
              FALSE               FALSE
                 NA                  NA
                 NA                  NA
                 NA                  NA
              FALSE               FALSE
              FALSE               FALSE
              FALSE               FALSE
               TRUE               FALSE
              FALSE               FALSE
              FALSE               FALSE
              FALSE               FALSE
                 NA                  NA
                 NA                  NA
                 NA                  NA
Finite donor-state factorial contrasts:
 policy year donor   delta_S   delta_I   Delta_S    Delta_I interaction reason
    raw 2015  2017  313251.4  14303.95 -8.213746 -10.403053    6.326963     NA
    raw 2017  2015 -313251.4 -14303.95 16.216159  10.672666   -8.769205     NA
  blend 2015  2017  301310.9  12622.62 -7.835502  -9.410372    5.964069     NA
  blend 2017  2015 -301310.9 -12622.62 15.742479   9.939037   -8.567231     NA

Initial-I timing effects are conditional fitted-model responses, not observed-history causal attribution.

The selected M2 endpoint is retained. Prior Phase 26 adequacy failures remain in docs/phase26_interpretation.md.
