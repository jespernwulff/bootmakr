* Runs behind "Why bootstrap?": the ITCV, analytic sensemakr and bootmakr on
* the simulated example data (firms.dta, copied into this folder by
* _prelude.do from BOOTMAKR_STATA_DIR).
* Start from this folder (see _prelude.do).
do "_prelude.do" why-bootstrap
quietly cd "$pc_out/.."

* >>> data
use firms, clear
* <<<

* >>> ses
quietly regress y x c
lincom x

quietly regress y x c, vce(cluster firm)
lincom x
* <<<

* >>> itcv
quietly regress y x c
scalar t_conv = _b[x] / _se[x]
scalar df = e(df_r)
quietly regress y x c, vce(cluster firm)
scalar t_rob = _b[x] / _se[x]

* the ITCV: t-implied correlation minus the critical one (Frank, 2000)
scalar r_crit = invttail(df, .025) / sqrt(invttail(df, .025)^2 + df)
foreach v in conv rob {
    scalar r_impl_`v' = t_`v' / sqrt(t_`v'^2 + df)
    scalar itcv_`v' = (r_impl_`v' - r_crit) / (1 - abs(r_crit))
}
display "conventional SE: t-implied correlation = " %5.3f r_impl_conv ///
    "   ITCV = " %5.3f itcv_conv
display "clustered SE:    t-implied correlation = " %5.3f r_impl_rob ///
    "   ITCV = " %5.3f itcv_rob "   sqrt(ITCV) = " %5.3f sqrt(itcv_rob)
* <<<

* >>> impact
* partial correlations given c, by residualising on c
foreach v in x y q {
    quietly regress `v' c
    predict double `v'_c, resid
}
quietly correlate x_c y_c
scalar r_xy_c = r(rho)
quietly correlate x_c q_c
scalar r_xq_c = r(rho)
quietly correlate y_c q_c
scalar r_yq_c = r(rho)
display "partial correlation of x and y given c: " %5.3f r_xy_c
display "q with x given c: " %5.3f r_xq_c "   q with y given c: " %5.3f r_yq_c ///
    "   impact of q: " %5.3f r_xq_c * r_yq_c
* <<<

* >>> withq
regress y x c q, vce(cluster firm)
* <<<

* >>> analytic
sensemakr y x c, treat(x) benchmark(c) kd(1)
* <<<

tic
* >>> bootmakr
bootmakr y x c, treat(x) benchmark(c) cluster(firm) seed(123)
* <<<
keep_result bootmakr
toc bootmakr

tic
* >>> rows
quietly bootmakr y x c, treat(x) benchmark(c) seed(123)
display "observations resampled:  SE = " %6.4f r(se) ///
    "   95% CI [" %6.4f r(ci_lower) ", " %6.4f r(ci_upper) "]"
* <<<
keep_result rows
toc rows

file close res
file close tim
log close pc
exit, clear
