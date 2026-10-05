* Runs behind "Get started in Stata", the Stata side of the home page and the
* comparison in "Clustered and stratified data". The data of Chen, Chittoor
* and Vissa (2021) are downloaded from the authors' OSF repository into this
* folder (ceo_pay.dta is not kept under version control).
* Start from this folder (see _prelude.do).
do "_prelude.do" first-steps
quietly cd "$pc_out/.."

* >>> data
global osf "https://osf.io/download/stk7y/?view_only=885d1ed5498c4b3c8a1194882e468286"
copy "$osf" "ceo_pay.dta", replace
use "ceo_pay.dta", clear
* <<<

* >>> controls
global controls "ceo_tenure PA_nic3_med ceo_edu_dummy lg_sales firm_age"
global controls "$controls promoters_pct institutions_pct i.year i.nic_code_1digit"
* <<<

* >>> regress
quietly regress lg_ceopay owner_ceo $controls, vce(cluster co_code)
lincom owner_ceo
* <<<

tic
* >>> first
bootmakr lg_ceopay owner_ceo $controls, ///
    treat(owner_ceo) benchmark(firm_age) cluster(co_code) seed(123)
* <<<
keep_result first

* >>> stored
return list
* <<<
toc first

tic
* >>> sweep
bootmakr lg_ceopay owner_ceo $controls, ///
    treat(owner_ceo) benchmark(firm_age) kd(0.5 1 1.5 2 2.5 3 3.5 4) ///
    cluster(co_code) seed(123) plot
* <<<
keep_sweep sweep 0.5 1 1.5 2 2.5 3 3.5 4

* >>> sweep-stored
matrix list r(results)
matrix list r(benchmark_strength)
* <<<
toc sweep
quietly graph export "$pc_out/first-steps-sweep.svg", replace

tic
* >>> program
capture program drop my_sens
program define my_sens, eclass
    sensemakr lg_ceopay owner_ceo $controls, ///
        treat(owner_ceo) benchmark(firm_age) kd(1) suppress
end

bootmakr, treat(owner_ceo) program(my_sens) cluster(co_code) seed(123)
* <<<
keep_result program
toc program

tic
* >>> compare
* firm-years resampled
quietly bootmakr lg_ceopay owner_ceo $controls, ///
    treat(owner_ceo) benchmark(firm_age) seed(123)
display "firm-years:  se = " %6.4f r(se) "   95% CI [" %6.4f r(ci_lower) ", " %6.4f r(ci_upper) "]   p = " %5.3f r(p)

* firms resampled
quietly bootmakr lg_ceopay owner_ceo $controls, ///
    treat(owner_ceo) benchmark(firm_age) cluster(co_code) seed(123)
display "firms:       se = " %6.4f r(se) "   95% CI [" %6.4f r(ci_lower) ", " %6.4f r(ci_upper) "]   p = " %5.3f r(p)
* <<<
toc compare

file close res
file close tim
log close pc
exit, clear
