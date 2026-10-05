* Worked example (Chen, Chittoor & Vissa, 2021), part 1 of 3: the data, the
* two sets of standard errors, the headline analysis (kd = 1), an omitted
* variable half as strong (kd = 0.5), and the seed check.
* Start from this folder; needs BOOTMAKR_CHEN_DIR (see _prelude.do).
do "_prelude.do" ceo-pay-1

* >>> setup
use "CEO pay project - Data for sharing FINAL.dta", clear

global controls "ceo_tenure PA_nic3_med ceo_edu_dummy lg_sales firm_age promoters_pct institutions_pct i.year i.nic_code_1digit"
* <<<

* >>> model2
quietly regress lg_ceopay owner_ceo $controls
lincom owner_ceo

quietly regress lg_ceopay owner_ceo $controls, vce(cluster co_code)
lincom owner_ceo
* <<<

tic
* >>> kd1
bootmakr lg_ceopay owner_ceo $controls, ///
    treat(owner_ceo) benchmark(ceo_tenure) cluster(co_code) ///
    reps(10000) seed(912323)
* <<<
keep_result kd1

* >>> kd1-stored
return list
* <<<
toc kd1

tic
* >>> kd05
bootmakr lg_ceopay owner_ceo $controls, ///
    treat(owner_ceo) benchmark(ceo_tenure) kd(0.5) cluster(co_code) ///
    reps(10000) seed(912323)
* <<<
keep_result kd05
toc kd05

tic
* >>> seeds
foreach s in 101 202 303 404 {
    quietly bootmakr lg_ceopay owner_ceo $controls, ///
        treat(owner_ceo) benchmark(ceo_tenure) cluster(co_code) seed(`s')
    display "seed `s':  p = " %5.3f r(p) "   95% CI [" %6.3f r(ci_lower) ", " %5.3f r(ci_upper) "]"
}
* <<<
toc seeds4x1000

file close res
file close tim
log close pc
exit, clear
