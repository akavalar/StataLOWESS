*! version 3.2.5  18jul2002
*! Fixed by Andraz (@akavalar), 13feb2017, no rights reserved

program define ksm_correct, sort
	version 7
	syntax varlist(min=2 max=2 numeric) [if] [in] [, /*
		*/ Line Weight noGraph GENerate(string) BWidth(real 0) /*
		*/ BY(varlist) Adjust LOWess LOgit T1title(string) /*
		*/ Symbol(string) Sort Connect(string) * ]
	local gen "`generate'"
	local generat

	if "`gen'"~="" {
		confirm new var `gen'
	}

	local tmpwgt "`weight'"
	local weight 
	marksample touse
	local weight "`tmpwgt'"
	local tmpwgt
	
	if "`lowess'"!="" {
		local weight "weight"
		local line "line"
	}

	local sort "sort"	/* Force -sort- */

	tokenize `varlist'
	local Y `1'
	local X `2'

	tempvar grp
	if "`by'"~="" {
		sort `touse' `by' 
		qui by `touse' `by': gen long `grp'=1 if _n==1 & `touse'
		qui replace `grp'=sum(`grp') if `touse'
		local n = `grp'[_N]
		local byopt "by(`by')"
					/* check more than one groups */
		if `n'<2 {
			di in red "by() variable takes on only one value"
			exit 198
		}
	}
	else {
		local n=1
		qui gen `grp'=1
	}
/*
	Set up complete.
*/
	local bw = `bwidth'
	if `bw' <= 0 | `bw' >= 1 { local bw .8 }

	local smlist
	local lsymbol
	local lconnect
	local i=1

	tempvar subuse
	qui gen byte `subuse' = .

	while `i'<=`n' {
		qui replace `subuse'=`touse' & `grp'==`i'
		tempvar Y1 X1 smooth`i'
		quietly {
			gen `Y1'=`Y' if `subuse'
			gen `X1'=`X' if `subuse'
			summ `X1'
			local cnt = r(N)
			if `cnt'<3 { noisily error 2001 }
			if r(max)-r(min)<1e-30 {
				di in red "Range of `2' is too small"
				exit 499
			}
			_crcslbl `Y1' `1'
			_crcslbl `X1' `2'
/*
	smopt=0 (unwtd mean), 1(unwtd line), 2(wtd mean), 3(wtd line).
*/
			local smopt = ("`line'"=="line") + cond("`weight'"=="weight",2,0)
			Ksm `Y1' `X1' `cnt' `bw' `smooth`i'' `smopt' `subuse'
/*
	Adjust smooth so that mean of Y-values equals mean of smoothed values.
*/
			if "`adjust'" == "adjust" {
				summ `smooth`i''
				local mean = r(mean)
				summ `Y1'
				replace `smooth`i'' = `smooth`i''* /*
					*/ r(mean)/`mean' if `subuse'
			}
			if "`logit'" == "logit" {
				local adj = 1/`cnt'
				local small = 0.0001
				replace `smooth`i'' = `adj' if `smooth`i'' /*
					*/ <`small' & `subuse'
				replace `smooth`i'' = 1-`adj'  /*
					*/ if `smooth`i''>(1-`small') & `subuse'
				replace `smooth`i'' = ln(`smooth`i''/(1- /*
					*/ `smooth`i'')) if `subuse'
			}
		}

		local smlist `smlist' `smooth`i''
		local lsymbol "`lsymbol'i"
		local lconnect "`lconnect'l"
		
		capture drop `Y1' `X1'
		local i=`i'+1
	}

/*
	Graph (if required)
*/
	if "`graph'" == "" { 
		if `"`t1title'"' ==""{
			if "`lowess'"!="" {
				local t1title "Lowess smoother"
			}
			else if "`line'" == "line" {
				local t1title "Running line smoother"
			}
			else 	local t1title "Running mean smoother"
			local t1title `"`t1title', bandwidth = `bw'"'
		}
		if "`logit'"=="" {
			tempvar sm
			qui egen `sm'=rmax(`smlist')
			if "`by'"~="" {sort `by'}

			local symbol = cond("`symbol'"=="", /*
				*/ "oi", "`symbol'")
			local connect = cond("`connect'"=="", /*
				*/ ".l", "`connect'")
			gr `Y' `sm' `X' if `touse', `options' /*
				*/ t1(`"`t1title'"') s(`symbol') /*
				*/ c(`connect') `sort' `byopt'
		}
		else {
			tempvar sm
			qui egen `sm'=rmax(`smlist')
			if "`by'"~="" {sort `by'}

			local symbol = cond("`symbol'"=="", /*
				*/ "i", "`symbol'")
			local connect = cond("`connect'"=="", /*
				*/ "l", "`connect'")
			gr `sm' `X' if `touse', `options' /*
				*/ t1(`"`t1title'"') s(`symbol') /*
				*/ c(`connect') `sort' `byopt'
		}
	}

	if "`gen'"~="" {
		qui egen `gen'=rmax(`smlist')
	}

end


program define Ksm
	args Y X count bwidth Ys smopt touse
	tempvar W
	capture drop `Ys'
/*
	Args: 1=Y, 2=X, 3=count, 4=bandwidth, 5=smoothed Y,
	6=0 (unwtd mean), 1(unwtd line), 2(wtd mean), 3(wtd line).
	NOTE: bandwidth is expressed as a fraction of the sample size,
	not a number of sample values.
*/

	tempvar revuse 
	gen byte `revuse' = cond(`touse',1,2)
	sort `revuse' `X'
	drop `revuse'
	/* touse sample is now 1/`count' */

/************************************ THIS IS NON-STANDARD, AVOID (ANDRAZ)
	Using symmetric neighbourhoods. 
	`k' is half-bandwidth. Initialise using points 1,...,k+1.
	Delta is distance between current Y and furthest neighbour in the
	interval.
*********************************** THIS IS NON-STANDARD, AVOID (ANDRAZ)
	
	`W' is tricube weight function within the interval.
*/
	

	gen double `W' = 0 in 1/`count'
	gen double `Ys' = .
	local wtd = (`smopt'>1)
	local mean = (`smopt'==0 | `smopt'==2)

	local ns = max(min(floor(`bwidth'*`count'),`count'),2) //size of neighborhood
	local i 0
	local nleft = 1 //starting point (left)
	local nright = `ns' //starting point (right)
	
	while `i'<`count' {
		local i = `i'+1
		local x0 = `X'[`i']
		
		* 3 possible cases:
		* 	- left corner: there are more points on the right of the point than there are on the left (asymmetric neighborhood); incrementing the point does not move the neighborhood
		* 	- middle case: there are equal number of points on both sides; incrementing the point moves the neighborhood to the right by one point
		* 	- right corner: there are more points on the left of the point than there are on the right (asymmetric neighborhood); incrementing the point does not move the neighborhood
		* Note that the size of the neighborhood is fixed at all times!
		
		local d1 = `x0' - `X'[`nleft']
		local d2 = `X'[`nright'+1] - `x0'
		
		//if incrementing the point means the left part of the neighborhood is bigger than the right one (d1 > d2), the neighborhood needs to be moved (middle case)
		//however, if `nright'==`count', we've reached the right corner and we cannot move the neighborhood anymore
		//if d1 <= d2, then do nothing as we're in the left corner and shouldn't move the neighborhood just yet
		if (`d1'>`d2') & (`nright'<`count') {
			local nleft = `nleft' + 1
			local nright = `nright' + 1
			}

		local h = max(`x0' - `X'[`nleft'], `X'[`nright'] - `x0')
		local h9 = 0.999*`h'
		local h1 = 0.001*`h'
		tempvar r
		quiet gen `r' = abs(`X'-`x0')
		
		if `wtd' {
			replace `W' = (1-(`r'/`h')^3)^3 if (`r'>`h1' & `r'<=`h9') in `nleft'/`nright' //all but the most distant value
			replace `W' = 1 if `r'<=`h1' in `nleft'/`nright' //weight=1 for the Y value corresponding to X[i]
			drop `r'
		}
		if `mean' {
			summ `Y' [aw=`W'] in `nleft'/`nright'
			replace `Ys' = r(mean) in `i'
		}
		else {
			reg `Y' `X' [aw=`W'] in `nleft'/`nright'
			replace `Ys' = _b[_cons]+_b[`X']*`x0' in `i'
		}
	}
end
