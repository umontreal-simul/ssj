Ordinary MC:
Asian option under GBM, for testing
REPORT on Tally stat. collector ==> value of Asian option
    num. obs.      min          max        average     variance    standard dev.
     100000        0.0000      279.2806       13.1572      526.0081       22.9349
  95.0% conf. interval for the mean (Student approx.): (    13.0151,    13.2994 )
Variance per run:   526.0081
Total CPU time:      0:0:0.9

------------------------

RQMC experiment:
Asian option under GBM, for testing
QMC point set: Sobol sequence:
DigitalNetBase2: , Num cols = 16, Num rows = 31, outDigits = 31
Randomization: Left matrix scramble + random digital shift
REPORT on Tally stat. collector ==> RQMC averages for Asian option under GBM
    num. obs.      min          max        average     variance    standard dev.
        50       13.109       13.136       13.122       2.8E-5       5.3E-3
  95.0% conf. interval for the mean (Student approx.): (    13.120,    13.123 )
Total CPU time:      0:0:1.89

----------------------------------------------------

Experiment for MC/RQMC comparison:
Asian option under GBM, for testing
REPORT on Tally stat. collector ==> Statistics with MC
    num. obs.      min          max        average     variance    standard dev.
     100000        0.0000      307.7575       13.0971      512.1530       22.6308
  95.0% conf. interval for the mean (Student approx.): (    12.9568,    13.2374 )
Variance per run:   512.1530
Total CPU time:      0:0:0.8

Asian option under GBM, for testing
QMC point set: Sobol sequence:
DigitalNetBase2: , Num cols = 16, Num rows = 31, outDigits = 31
Randomization: Left matrix scramble + random digital shift
REPORT on Tally stat. collector ==> Statistics on RQMC averages
    num. obs.      min          max        average     variance    standard dev.
        50       13.114       13.133       13.122       1.8E-5       4.3E-3
  95.0% conf. interval for the mean (Student approx.): (    13.121,    13.123 )
Total CPU time:      0:0:1.89

MC Variance per run:    512.15303
RQMC Variance per run:    1.19452
Variance ratio:        428.752
Efficiency ratio:      580.552
---------------------------------------------------------------
