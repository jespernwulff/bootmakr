* Runs behind "Get started in Stata": the simulated example data that come
* with the command (firms.dta, copied into this folder by _prelude.do from
* BOOTMAKR_STATA_DIR; not kept under version control).
* Start from this folder (see _prelude.do).
do "_prelude.do" getting-started
quietly cd "$pc_out/.."

* >>> data
use firms, clear
describe
* <<<

* >>> naive
regress y x c, vce(cluster firm)
* <<<

tic
* >>> first
bootmakr y x c, treat(x) benchmark(c) cluster(firm) seed(123)
* <<<
keep_result first

* >>> stored
return list
* <<<
toc first

tic
* >>> sweep
bootmakr y x c, treat(x) benchmark(c) kd(0.5 0.75 1 1.25 1.5) ///
    cluster(firm) seed(123) plot
* <<<
keep_sweep sweep 0.5 0.75 1 1.25 1.5

* >>> sweep-stored
matrix list r(results)
matrix list r(benchmark_strength)
* <<<
toc sweep
quietly graph export "$pc_out/getting-started-sweep.svg", replace

* >>> check
regress y x c q, vce(cluster firm)
* <<<

tic
* >>> program
capture program drop my_sens
program define my_sens, eclass
    sensemakr y x c, treat(x) benchmark(c) kd(1) suppress
end

bootmakr, treat(x) program(my_sens) cluster(firm) seed(123)
* <<<
keep_result program
toc program

file close res
file close tim
log close pc
exit, clear
