* Convert the help file of the Stata command to plain text for the website
* article "Stata help file". Start from this folder; BOOTMAKR_STATA_DIR must
* point to the folder that holds bootmakr.sthlp.
clear all
set more off
local adodir : environment BOOTMAKR_STATA_DIR
if "`adodir'" == "" {
    display as error "set BOOTMAKR_STATA_DIR to the folder with bootmakr.sthlp"
    exit 198
}
capture mkdir "output"
translator set smcl2txt linesize 92
quietly translate "`adodir'/bootmakr.sthlp" "output/bootmakr-help.txt", ///
    translator(smcl2txt) replace
exit, clear
