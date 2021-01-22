use experiment_hiD1.20_collapse082019.dta, clear
gen durables_diff = nextdurables - durables
twoway (line durables_diff age if marginal==0) (line durables_diff age if marginal==1) if income_bin >= 4 & income_bin <= 12, yline(0, lc(gs2)) legend(label(2 "Marginal") label(1 "Inframarginal")) by(income_bin) name(withmin_agg, replace)
use experiment_monetary_collapse082019.dta, clear
gen durables_diff = nextdurables - durables
twoway (line durables_diff age if marginal==0) (line durables_diff age if marginal==1) if income_bin >= 6 & income_bin <= 12, yline(0, lc(gs2)) legend(label(2 "Marginal") label(1 "Inframarginal")) by(income_bin) name(nomin_agg, replace)
use experiment_hiD1.20_collapse082019.dta, clear
gen durables_diff = nextdurables - durables
twoway (line durables_diff age if marginal==0) (line durables_diff age if marginal==1) if income_bin >= 4 & income_bin <= 12, yline(0, lc(gs2)) legend(label(2 "Marginal") label(1 "Inframarginal")) by(income_bin) name(withmin_agg, replace)
twoway (line nextdurables age if marginal==0) (line nextdurables age if marginal==1) if income_bin >= 4 & income_bin <= 12, yline(0, lc(gs2)) legend(label(2 "Marginal") label(1 "Inframarginal")) by(income_bin) name(withmin_levels, replace)
twoway (line nextdurables age if marginal==0) (line nextdurables age if marginal==1) if income_bin >= 4 & income_bin <= 12, yline(1.2, lc(gs2)) legend(label(2 "Marginal") label(1 "Inframarginal")) by(income_bin) name(withmin_levels, replace)
use experiment_monetary_collapse082019.dta, clear
twoway (line nextdurables age if marginal==0) (line nextdurables age if marginal==1) if income_bin >= 4 & income_bin <= 12, yline(0.25, lc(gs2)) legend(label(2 "Marginal") label(1 "Inframarginal")) by(income_bin) name(nomin_levels, replace)
twoway (line nextdurables age if marginal==0) (line nextdurables age if marginal==1) if income_bin >= 6 & income_bin <= 12, yline(0.25, lc(gs2)) legend(label(2 "Marginal") label(1 "Inframarginal")) by(income_bin) name(nomin_levels, replace)
twoway (line nextdurables age if marginal==0) (line nextdurables age if marginal==1) if income_bin >= 6 & income_bin <= 12, yline(1.2, lc(gs2)) legend(label(2 "Marginal") label(1 "Inframarginal")) by(income_bin) name(nomin_levels, replace)
tab income_bin marginal [iw=id]
tab income_bin marginal [iw=id], column
use experiment_hiD1.20_collapse082019.dta, clear
tab income_bin marginal [iw=id], column
