* Runs behind the Stata side of the home page, the data block of the worked
* example and the comparison in "Clustered and stratified data". The data of
* Chen, Chittoor and Vissa (2021) are downloaded from the authors' OSF
* repository into this folder (ceo_pay.dta is not kept under version
* control). Start from this folder (see _prelude.do).
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
