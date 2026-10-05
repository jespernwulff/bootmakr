* Shared set-up for the Stata runs behind the website articles. Not shown on
* the site. Each run file calls it as
*     do "_prelude.do" <name of the run>
* and is started from this folder (precompute/stata).
*
* Environment variables:
*   BOOTMAKR_STATA_DIR  folder holding bootmakr.ado, if the version to document
*                       is not the installed one
*   BOOTMAKR_CHEN_DIR   folder holding the Chen et al. (2021) data file; needed
*                       by the ceo-pay runs only
args run

clear all
set more off
set graphics off
set linesize 120

* Everything a run produces goes to output/: the log the articles quote from
* (<run>.txt), key results, timings and graphs. Stata's own batch log,
* <run>.log in this folder, starts with the licence banner and is not kept.
capture mkdir "output"
global pc_out "`c(pwd)'/output"

local adodir : environment BOOTMAKR_STATA_DIR
if "`adodir'" != "" {
    adopath ++ "`adodir'"
}
capture program drop bootmakr

* Key results of each bootmakr call go to <run>-results.txt so that the
* article text can quote them without retyping.
capture file close res
file open res using "$pc_out/`run'-results.txt", write text replace
file write res "label" _tab "kd" _tab "estimate" _tab "se" _tab "ci_lower" _tab "ci_upper" _tab "p" ///
    _tab "N" _tab "N_reps" _tab "N_successful" _n

capture program drop keep_result
program define keep_result
    * usage: keep_result label   (after a bootmakr call with a single kd)
    args label
    file write res "`label'" _tab "." _tab %12.0g (r(estimate)) _tab %12.0g (r(se)) ///
        _tab %12.0g (r(ci_lower)) _tab %12.0g (r(ci_upper)) _tab %12.0g (r(p)) ///
        _tab %12.0g (r(N)) _tab %12.0g (r(N_reps)) _tab %12.0g (r(N_successful)) _n
end

capture program drop keep_sweep
program define keep_sweep
    * usage: keep_sweep label kd1 kd2 ...   (after a bootmakr call with several kd)
    gettoken label 0 : 0
    tempname R
    matrix `R' = r(results)
    local N = r(N)
    local N_reps = r(N_reps)
    local N_ok = r(N_successful)
    local i = 0
    foreach k of local 0 {
        local ++i
        file write res "`label'" _tab "`k'" _tab %12.0g (`R'[`i', 1]) _tab %12.0g (`R'[`i', 2]) ///
            _tab %12.0g (`R'[`i', 3]) _tab %12.0g (`R'[`i', 4]) _tab %12.0g (`R'[`i', 5]) ///
            _tab %12.0g (`N') _tab %12.0g (`N_reps') _tab %12.0g (`N_ok') _n
    }
end


* Timing: tic ... toc label  (writes seconds to <run>-timings.txt). toc runs
* -timer list-, which replaces r(): call keep_result or keep_sweep first.
capture file close tim
file open tim using "$pc_out/`run'-timings.txt", write text replace
file write tim "label" _tab "seconds" _n
capture program drop tic
program define tic
    timer clear 99
    timer on 99
end
capture program drop toc
program define toc
    args label
    timer off 99
    quietly timer list 99
    file write tim "`label'" _tab %9.1f (r(t99)) _n
end

* Convergence summary of a bootmakr call with converge().
capture program drop keep_conv
program define keep_conv
    args label
    capture file close cnv
    quietly file open cnv using "$pc_out/`label'-convergence.txt", write text replace
    file write cnv "se_mean" _tab "se_range" _tab "se_cv" _tab "se_range_high" _tab "se_cv_high" ///
        _tab "p_mean" _tab "p_range" _tab "p_cv" _tab "p_range_high" _tab "p_cv_high" _tab "threshold" _n
    file write cnv %12.0g (r(conv_se_mean)) _tab %12.0g (r(conv_se_range)) _tab %12.0g (r(conv_se_cv)) ///
        _tab %12.0g (r(conv_se_range_high)) _tab %12.0g (r(conv_se_cv_high)) ///
        _tab %12.0g (r(conv_p_mean)) _tab %12.0g (r(conv_p_range)) _tab %12.0g (r(conv_p_cv)) ///
        _tab %12.0g (r(conv_p_range_high)) _tab %12.0g (r(conv_p_cv_high)) _tab %12.0g (r(conv_thresh_reps)) _n
    file close cnv
end

* The run files read the Chen et al. (2021) data by file name, as a reader
* would, so they are started in the folder that holds it.
local chendir : environment BOOTMAKR_CHEN_DIR
if "`chendir'" != "" {
    quietly cd "`chendir'"
}

* The log that the articles quote from. It is opened last: after Stata's
* start-up banner (no licence details) and after the set-up above (no local
* paths).
capture log close pc
log using "$pc_out/`run'.txt", text replace nomsg name(pc)
