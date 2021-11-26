
//local drive Z:/
local drive /Users/Michael/Documents
cd "`drive'/merlin"
adopath ++ "`drive'/merlin"
adopath ++ "`drive'/merlin/merlin"
clear all

do ./build/buildmlib.do
mata mata clear

local drive /Users/Michael/Documents/stmixed

cd "`drive'"
adopath ++ "`drive'/stmixed"
clear all

set seed 130931

clear
set obs 4
gen id1 = _n
expand 300
gen trt = runiform()>0.5
bys id1 : gen u1 = rnormal(0,1)
bys id1 : replace u1 = u1[1]

survsim stime dead , dist(weib) lambda(0.1) gamma(1.2) cov(trt -0.5 u1 1) maxt(5)

stset stime, f(dead)

// mestreg trt || id1: , dist(weib) nohr //intmethod(gh) intpoints(15)
// merlin (_t trt M1[id1]@1, family(weibull, failure(_d))), 

stmixed trt || id1: , distribution(rp) df(3) 
stmixed trt || id1: , distribution(rcs) df(3) 
stmixed trt || id1: , distribution(pwe) knots(1 2 4)

foreach dist in exp weib gomp logn logl {
	stmixed trt || id1: , distribution(`dist') 
}
