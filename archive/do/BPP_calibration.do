/*******************************************************************************
    Title: BPP_calibration.do
    Purpose: Clean SCF data for model moment calibration.
    Author: TC, DB
    Date: 2016-03-31

*******************************************************************************/

/* There is very little code in this file.
   The way the procedure goes is that we use BPP's replication files,
   specifically mindist_AER.do, all the way up to when log income is residualized.
   Our definition of income is household earnings and before taxes
   (add them in the model), so the code is 500000.
   and remove all "year of birth" and "education" dummies from the dummies regression.

   The leftover residuals, along with the age of the interviewee, are saved
   into dataBPP.dta. Now we fit a parametric function of age on them.

*/

capture program drop make_BPP
program define make_BPP

    use "$dboxdatadir/dataBPP.dta", clear
    keep if age >= 22 & age < 60 // Retirement happens ON YEAR 60

    binscatter uy age, n(40) title("BEFORE taxes")
    sleep 3000
    forvalues i=2/4 {
        gen age`i' = age^`i'
    }
    * Cubic fit. note the if condition: see later note
    reg uy age age2-age3  if age >= 25
    predict logy_model_wo, xb

    svyset [pw=weight]
    svy: mean uy
    matrix b = e(b)
    scalar means = b[1,1]
    display exp(means)

    collapse (first) logy_model (mean) uy (p75) uy_upper=uy (p25) uy_lower=uy, ///
        by(age)
    * The cubic fit overestimates the earliest years, so for now just linearly
    * extrapolate over them? In future may just start at age 22/23
    replace logy_model = uy if inlist(age, 23)
    replace logy_model = (uy[2] + uy[4])/2 if age == 24
    replace logy_model = logy_model[2] - (logy_model[4]-logy_model[3]) if age == 22

end

make_BPP
* Diagnostic: how well does the polynomial fit?
twoway scatter uy logy_model_wo age, title("BEFORE taxes")
sleep 3000

replace logy_model_wo = logy_model_wo - means
format logy_model_wo %13.10f
keep logy_model_wo
export delimited $initdir/ageearnings.txt, novar dataf replace

