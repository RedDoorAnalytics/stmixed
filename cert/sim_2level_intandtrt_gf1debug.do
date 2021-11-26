//2 level, random intercept and trt, check score on true model

local drive Z
cd "`drive':\stdev\stmixed"
adopath ++ ".\stmixed"
clear all

do ./build/buildmlib.do

set seed 7254

clear
set obs 500 
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
stmixed trt || id1: trt, d(w) nonadapt gh(15) mlmethod(gf1) cov(exch)


stmixed trt || id1: trt, d(w) nonadapt gh(15) mlmethod(gf1debug) cov(indep)
stmixed trt || id1: trt, d(w) nonadapt gh(15) mlmethod(gf1debug) cov(unstr)
stmixed trt || id1: trt, d(w) nonadapt gh(15) mlmethod(gf1debug) cov(iden)
