* Worked example (Chen, Chittoor & Vissa, 2021), part 2 of 3: the breakdown
* point (kd sweep with firm age as benchmark) and a grouped benchmark.
* Start from this folder; needs BOOTMAKR_CHEN_DIR (see _prelude.do).
do "_prelude.do" ceo-pay-2

use "CEO pay project - Data for sharing FINAL.dta", clear
global controls "ceo_tenure PA_nic3_med ceo_edu_dummy lg_sales firm_age promoters_pct institutions_pct i.year i.nic_code_1digit"

tic
* >>> sweep
bootmakr lg_ceopay owner_ceo $controls, ///
    treat(owner_ceo) benchmark(firm_age) kd(0.5 1 1.5 2 2.5 3 3.5 4) ///
    cluster(co_code) reps(10000) seed(912323) plot
* <<<
keep_sweep sweep 0.5 1 1.5 2 2.5 3 3.5 4
toc sweep
quietly graph export "$pc_out/ceo-pay-sweep.svg", replace

tic
* >>> grouped
bootmakr lg_ceopay owner_ceo $controls, ///
    treat(owner_ceo) gbenchmark(ceo_tenure ceo_edu_dummy) kd(0.1 0.2 0.3 0.4 0.5) ///
    cluster(co_code) reps(10000) seed(912323) plot
* <<<
keep_sweep grouped 0.1 0.2 0.3 0.4 0.5
toc grouped
quietly graph export "$pc_out/ceo-pay-grouped.svg", replace

file close res
file close tim
log close pc
exit, clear
