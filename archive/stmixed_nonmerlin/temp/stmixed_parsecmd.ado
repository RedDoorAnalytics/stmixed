*! version 1.0.0 ?????2015

/*
History
version 1.0.0 - adapted from _xtmixed_parsecmd.ado
*/

program stmixed_parsecmd, sclass
	version 13.1
	gettoken i 0 : 0

	// First we strip off if and in - they have already been processed
	// global weights too
	
	syntax [anything] [if] [in] [, * ]
	local 0 `anything', `options'

	local myopts noConstant COVariance(name) COLlinear 

	gettoken lev rest : 0 , parse(":")

	gettoken colon rest : rest , parse(":")
	if `"`colon'"' == `":"' {
		local lev = trim("`lev'")
		if "`lev'" != "_all" {
			unab lev : `lev'
			confirm variable `lev'
		}
		else {
			di as error "_all not allowed" 
			exit 198
		}
		local 0 `"`rest'"'
		capture syntax [varlist(numeric default=none)] [,*]
		if _rc == 101 {
			di as err "{p 0 6 2}must use R. when specifying "
			di as err "factor variables in random-effects "
			di as err "equations{p_end}"
			exit 198
		}
	}
	else {
		if (`i' > 0)  {
			di as err "{p 0 6 2}invalid random-effects "
			di as err "specification; perhaps you omitted the "
			di as err "colon after the level "
			di as err "variable{p_end}"
			exit 198
		}
		else	      local lev
	}

	capture syntax [varlist(numeric default=none)]		///
		       [ , `myopts' ]

	if _rc {						// r.varname or anything invalid
	    if `i' > 0 {
		    di as error `"`anything' invalid level specification"'
		    exit 198
		}
	}

	ParseCov cov : , `covariance'

	if "`cov'" != "identity" & `:list sizeof varlist' + ("`constant'"=="") == 1 {
	    if `"`covariance'"' != `""' {
			if "`lev'"!="" {
				local name "in {inp:`lev'} equation"
			}
			local msg di as text "{p 0 6} Note: single-variable random-"
			local msg `msg' "effects specification `name'; covariance "
			local msg `msg' "structure set to {inp:identity}{p_end}"		
			sreturn local msg `msg'
	    }
	    local cov identity
	}

	if `i' == 0  {					// parsing fixed level
	    if `"`lev'"' != "" {
		di as err "{p 0 4 2}level specification not allowed "
		di as err "for fixed equation{p_end}"
		exit 198
	    }

	    local cov identity

	    if "`covariance'" != "" {
			di as error						///
			"option covariance() not allowed with fixed-effects equation"
			exit 198
	    }
	    if "`collinear'" != "" {
	    	di as error						///
			"option collinear not allowed with fixed-effects equation"
			exit 198
	    }

	}

	c_local levnm_`i' `lev'
	c_local cov_`i'      `cov'
	c_local constant_`i' "`constant'"
	c_local collin_`i' "`collinear'"
	c_local varnames_`i'  `varlist'
	
end

program ParseCov
	gettoken lmac   0 : 0
	gettoken colon  0 : 0

	syntax [ , INDependent EXchangeable IDentity UNstructured * ]

	if `"`options'"' != "" {
		di as error `"`options' invalid covariance structure"'
		exit 198
	}

	local res `independent' `exchangeable' `identity' `unstructured'

	if ("`res'" == "") local res "independent"

	c_local `lmac' `res'
end

exit
