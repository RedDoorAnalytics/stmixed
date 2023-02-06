
local drive /Users/Michael/Documents/reddooranalytics/products
cd "`drive'/merlin"
adopath ++ "`drive'/merlin"
adopath ++ "`drive'/merlin/merlin"
clear all

do ./build/buildmlib.do
mata mata clear

cd "`drive'/stmixed"
adopath ++ "`drive'/stmixed"
adopath ++ "`drive'/stmixed/stmixed"
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

stmixed  || id1: , distribution(rp) df(3) 
stmixed  || id1: , distribution(rcs) df(3) 
stmixed trt || id1: , distribution(rp) df(3) 
stmixed trt || id1: , distribution(rcs) df(3) 
// stmixed trt || id1: , distribution(pwe) knots(1 2 4)

foreach dist in exp weib gomp logn logl {
        stmixed  || id1: , distribution(`dist') 
	stmixed trt || id1: , distribution(`dist') 
}
