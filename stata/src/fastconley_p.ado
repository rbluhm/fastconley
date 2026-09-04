*! version 0.2.0 04sep2026  predict after fastconley IV/2SLS
program define fastconley_p, rclass
	version 14.1
	syntax newvarname [if] [in] [, XB Residuals]
	opts_exclusive "`xb' `residuals'"
	if ("`xb'" == "" & "`residuals'" == "") {
		di as text "(option xb assumed; fitted values)"
		loc xb xb
	}
	if ("`residuals'" != "") {
		if ("`e(resid)'" == "") {
			di as error "residual prediction requires residuals(newvar) on the fastconley estimation command"
			exit 198
		}
		confirm numeric variable `e(resid)', exact
		gen double `varlist' = `e(resid)' `if' `in'
		label variable `varlist' "Residuals"
		exit
	}
	_predict double `varlist' `if' `in', xb
	label variable `varlist' "Linear prediction"
end
