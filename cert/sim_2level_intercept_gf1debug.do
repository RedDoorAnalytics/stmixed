//2 level, random intercept, check score on true model

local drive Z

cd "`drive':\stdev\stmixed"
adopath ++ ".\stmixed"
clear all

do ./build/buildmlib.do

set seed 7254

clear
set obs 500
gen id1 = _n
expand 20
bys id1: gen id2 = _n
gen trt = runiform()>0.5
bys id1 (id2) : gen u1 = rnormal() if _n==1
bys id1 (id2) : replace u1 = u1[1]

survsim stime dead , dist(weib) lambda(0.1) gamma(1.2) cov(trt -0.5 u1 1) maxt(5)
stset stime, f(dead)

stmixed trt || id1:, dist(weib) gh(15) nonadapt mlmethod(gf1debug)
