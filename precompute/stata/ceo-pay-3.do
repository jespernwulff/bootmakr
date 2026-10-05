* Worked example (Chen, Chittoor & Vissa, 2021), part 3 of 3: the convergence
* check at kd = 1 with 40,000 replications.
* Start from this folder; needs BOOTMAKR_CHEN_DIR (see _prelude.do).
do "_prelude.do" ceo-pay-3

use "CEO pay project - Data for sharing FINAL.dta", clear
global controls "ceo_tenure PA_nic3_med ceo_edu_dummy lg_sales firm_age promoters_pct institutions_pct i.year i.nic_code_1digit"

tic
* >>> converge
bootmakr lg_ceopay owner_ceo $controls, ///
    treat(owner_ceo) benchmark(ceo_tenure) cluster(co_code) ///
    reps(40000) seed(912323) ///
    converge(minreps(1000) stepsize(1000) threshold(10000))
* <<<
keep_conv ceo-pay-3
keep_result converge
toc converge
quietly graph export "$pc_out/ceo-pay-converge.svg", replace

file close res
file close tim
log close pc
exit, clear
