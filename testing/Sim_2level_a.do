
local drive /Users/Michael/Documents/reddooranalytics/products/stmixed

cd "`drive'"
adopath ++ "`drive'/stmixed"
clear all

set seed 7254

clear
set obs 1000 
gen id1 = _n
expand 10
bys id1: gen id2 = _n
gen trt = runiform()>0.5
mat cor = (1,0\0,1)
drawnorm u1 u2, means(0 0) sds(1 1) corr(cor)
//bys id1 (id2) : gen u1 = rnormal() if _n==1
bys id1 (id2) : replace u1 = u1[1]
//bys id1 (id2) : gen u2 = rnormal() if _n==1
bys id1 (id2) : replace u2 = u2[1]
gen trtui = (-0.5+u2) * trt

survsim stime dead , dist(weib) lambda(0.1) gamma(1.2) cov(trtui 1 u1 1) maxt(5)
stset stime, f(dead)

gen bhaz = 0

//mestreg trt || id1: trt, dist(weib) intpoints(7) 

// merlin (_t trt M1[id1]@1 trt#M2[id1]@1 , family(weibull, failure(_d))), 

stmixed trt || id1: trt, distribution(addrcs) df(3) intpoints(9) bhazard(bhaz) tvc(trt) dftvc(1) adaptopts(log)

