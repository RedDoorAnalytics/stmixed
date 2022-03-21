//build new version of merlin
//either for website or SSC
// --> run whole do file

local drive /Users/Michael/Documents/reddooranalytics/products/stmixed
cd `drive'

local sscbuild 		= 1

//=======================================================================================================================//

//build for SSC -> current version up is 2.2.1
if `sscbuild' {
	local sscversion 2_2_1
	cap mkdir ./ssc/version_`sscversion'
	local fdir `drive'/ssc/version_`sscversion'/
}


//=======================================================================================================================//

//pkg files
if `sscbuild' {
	copy ./ssc/stmixed_details.txt `fdir', replace
}
else {
	copy ./stmixed/stmixed.pkg `fdir', replace
	copy ./stmixed/stata.toc `fdir', replace
}
	
//=======================================================================================================================//

//stmixed
copy ./stmixed/stmixed.ado `fdir', replace

//help files
copy ./stmixed/stmixed.sthlp `fdir', replace
copy ./stmixed/stmixed_postestimation.sthlp `fdir', replace

