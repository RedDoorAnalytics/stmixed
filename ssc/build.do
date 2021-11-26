//build new version of merlin
//either for website or SSC
// --> run whole do file


//local drive Z:/
local drive /Users/Michael/Documents
cd `drive'/stmixed/

local sscbuild 		= 0

//=======================================================================================================================//

//build for SSC -> current version up is 2.0.3
if `sscbuild' {
	local sscversion 2_1_0
	cap mkdir ./ssc/version_`sscversion'
	local fdir /Users/Michael/Documents/stmixed/ssc/version_`sscversion'/
}
//build for website -> 2.2.0
else {
	local fdir /Users/Michael/Documents/website/static/code/stmixed/
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

