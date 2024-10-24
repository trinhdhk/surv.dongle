# Aalen-Johansen estimator survival plots for competing risks survival analysis

This dongle includes two functions for tidying and plotting AJ estimator survival curves for competing risks. 

1. `tidy_ajcomprsk`: Create tidied data tailored plotting AJ curves.
2. `gg_ajsurvplot` and `gg_ajsurvplot2`: create AJ cum-haz curve for the event of interest, without and with one competing event. The former supports both `survminer::ggsurvplot` and `ggsurvfit::ggsurvfit` as backend with minimal change to their parameters.

3. Trinh Dong, 2024
