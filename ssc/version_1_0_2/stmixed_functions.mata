*! version 0.2.0 13may2013 MJC

// stmixed Mata functions

//panelsum() added

/*
History
MJC 13may2013 version 0.2.0 - AFTs added
MJC 25jul2012 version 0.1.0
*/

version 12.1

local SS 	string scalar
local RS	real scalar
local NM	numeric matrix

mata:

//d0 core function
void stmixed_d0(	transmorphic scalar M,
					real scalar todo,
					real rowvector b,
					real scalar lnf,
					real rowvector g,
					real matrix H)
{
	pointer(struct stmixed_main scalar) scalar S
	S = &moptimize_util_userinfo(M,1)

	(*(*S).fillbetaP)(S,M,b)
	lnf = moptimize_util_sum(M,(*(*S).lnlP)(S))
	if (todo==0) return
}

//refine starting values with non-adaptive quadrature
void stmixed_d0_refine(	transmorphic scalar M,
						real scalar todo,
						real rowvector b,
						real scalar lnf,
						real rowvector g,
						real matrix H)
{
	pointer(struct stmixed_main scalar) scalar S
	S = &moptimize_util_userinfo(M,1)

	(*(*S).fillbetaP)(S,M,b)
	lnf = moptimize_util_sum(M,(*(*S).refinelnlP)(S))
	if (todo==0) return
}

// Exponential
void stmixed_e_bmat(pointer(struct stmixed_main scalar) scalar S, 
					scalar M,
					rowvector b)
{
	main_l_mat = exp(moptimize_util_xb(M,b,1))
	if ((*S).hascovs) S->l_mat = main_l_mat
	else S->l_mat = J((*S).Nobs,1,main_l_mat)
	stmixed_fillvcv(S,M,b,2)
}

// Weibull
void stmixed_w_bmat(pointer(struct stmixed_main scalar) scalar S, 
					scalar M,
					rowvector b)
{
	main_l_mat = exp(moptimize_util_xb(M,b,1))
	if ((*S).hascovs) S->l_mat = main_l_mat
	else S->l_mat = J((*S).Nobs,1,main_l_mat)
	S->g1 = exp(moptimize_util_xb(M,b,2))
	stmixed_fillvcv(S,M,b,3)
}

// Gompertz
void stmixed_g_bmat(pointer(struct stmixed_main scalar) scalar S, 
					scalar M,
					rowvector b)
{
	main_l_mat = exp(moptimize_util_xb(M,b,1))
	if ((*S).hascovs) S->l_mat = main_l_mat
	else S->l_mat = J((*S).Nobs,1,main_l_mat)
	S->g1 = moptimize_util_xb(M,b,2)
	stmixed_fillvcv(S,M,b,3)
}

// FPM
void stmixed_fpm_bmat(	pointer(struct stmixed_main scalar) scalar S, 
						scalar M,
						rowvector b)
{
	S->l_mat = moptimize_util_xb(M,b,1)
	S->dxb 	 = moptimize_util_xb(M,b,2)
	stmixed_fillvcv(S,M,b,3)
}

// AFTs
void stmixed_aft_bmat(pointer(struct stmixed_main scalar) scalar S, 
						scalar M,
						rowvector b)
{
	if ((*S).llog) main_l_mat = exp(-moptimize_util_xb(M,b,1))
	else main_l_mat = moptimize_util_xb(M,b,1)
	if ((*S).hascovs) S->l_mat = main_l_mat
	else S->l_mat = J((*S).Nobs,1,main_l_mat)
	S->g1 = exp(moptimize_util_xb(M,b,2))
	if ((*S).gam) {
		S->kapp = moptimize_util_xb(M,b,3)
		stmixed_fillvcv(S,M,b,4)
	}
	else stmixed_fillvcv(S,M,b,3)
}

// fill variance-covariance parameters
void stmixed_fillvcv(	pointer(struct stmixed_main scalar) scalar S, 
						scalar M,
						rowvector b,
						`RS' eq_ind)
{
	if ((*S).covind | (*S).covunstr) {
		for(i=1;i<=(*S).Nres;i++) {
			S->vcv[i,i] =  exp(moptimize_util_xb(M,b,eq_ind)):^2
			eq_ind++
		}
	}
	else {
		var_xb = exp(moptimize_util_xb(M,b,eq_ind)):^2
		for(i=1;i<=(*S).Nres;i++) S->vcv[i,i] = var_xb
		eq_ind++
	}
	if ((*S).covexch & (*S).Nres>1) {
		corr_xb = tanh(moptimize_util_xb(M,b,eq_ind))
		eq_ind = eq_ind :+ 1
		ind = 1
		while(ind<(*S).Nres){
			for(i=ind:+1;i<=(*S).Nres;i++){
				S->vcv[ind,i] = S->vcv[i,ind] = var_xb:*corr_xb
			}
			ind++
		}
	}		
	else if ((*S).covunstr) {
		ind = 1
		while(ind<(*S).Nres){
			for(i=ind:+1;i<=(*S).Nres;i++){
				corr_xb = tanh(moptimize_util_xb(M,b,eq_ind))
				S->vcv[ind,i] = S->vcv[i,ind] = sqrt((*S).vcv[ind,ind]):*sqrt((*S).vcv[i,i]):*corr_xb
				eq_ind++
			}
			ind++
		}		
	}
}

// Update adaptive quadrature locations and scales
void stmixed_update_quad(pointer(struct stmixed_main scalar) scalar SS)
{
	like_j = (*SS).detcholVCVs :* (*SS).jlnodes * (*SS).adaptbaseGHweights
	Nres = (*SS).Nres
	
	for(i=1;i<=(*SS).N;i++) {
		nodes_i = asarray((*SS).GHNarray,i)	
		test1 = (*SS).detcholVCVs[i,] :* (*SS).jlnodes[i,]:/like_j[i,] 
		newblups =  test1 * (nodes_i :* (*SS).adaptbaseGHweights')'
		vcv_new = (test1 * (stmixed_outerprod_by_col(nodes_i') :* (*SS).adaptbaseGHweights)) :- rowshape((newblups')*newblups,1)
		vcv_new = rowshape(vcv_new,Nres)
		stmixed_seblup_check(min(diagonal(vcv_new)))	//error check
		nodes = newblups' :+ cholesky(vcv_new) * (*SS).baseGHnodes
		asarray(SS->GHNarray,i,nodes)
		SS->detcholVCVs[i,] = sqrt(det(vcv_new))
		if (i==1) testmat = nodes'
		else testmat = testmat\nodes'
	}	
	SS->stackednodes = testmat
}	

// prolog
void stmixed_prolog(real rowvector b,scalar M,real scalar lnl) 
{
	pointer(struct stmixed_main scalar) scalar pSTJM		
	pragma unused lnl
	pSTJM = &moptimize_util_userinfo(M,1)
	
	atol = (*pSTJM).atol
	
	oldlnl = lnl
	if ((*pSTJM).iter==0) {
		(*(*pSTJM).fillbetaP)(pSTJM,M,b)
		oldlnl = moptimize_util_sum(M,(*(*pSTJM).lnlP)(pSTJM))
	}
	n = 0
	st_numscalar("mliter",n)
	st_numscalar("ll",oldlnl)		
	if ((*pSTJM).showadapt) stata(`"di as txt "-- Iteration " _col(14) mliter as txt ":" _col(19) "Adapted log likelihood = " as res %10.0g ll"')
	
	n = 1
	stmixed_update_quad(pSTJM)
	newlnl = moptimize_util_sum(M,(*(*pSTJM).lnlP)(pSTJM))
	st_numscalar("mliter",n)
	st_numscalar("ll",newlnl)		
	if ((*pSTJM).showadapt) stata(`"di as txt "-- Iteration " _col(14) mliter as txt ":" _col(19) "Adapted log likelihood = " as res %10.0g ll"')

	while (abs(oldlnl:-newlnl)>atol & n<(*pSTJM).reiterate) {
		n++
		swap(oldlnl,newlnl)
		stmixed_update_quad(pSTJM)
		newlnl = moptimize_util_sum(M,(*(*pSTJM).lnlP)(pSTJM))
		st_numscalar("mliter",n)
		st_numscalar("ll",newlnl)		
		if ((*pSTJM).showadapt) stata(`"di as txt "-- Iteration " _col(14) mliter as txt ":" _col(19) "Adapted log likelihood = " as res %10.0g ll"')
	}
	
	if (n==(*pSTJM).reiterate) {
		errprintf("Adaptive quadrature iterations have failed to converge within "+strofreal((*pSTJM).reiterate)+" iterations.\n")
		errprintf("Try reducing adaptit() or non-adaptive quadrature.\n")
		exit(1986)
	}
	
	pSTJM->iter = (*pSTJM).iter :+ 1
}

// multivariate normal probability density function
`NM' stmixed_dmvnorm(`NM' X, `NM' V)
{
	`NM' invsigma, Vvecs, Vvals, q2, denom
	
	symeigensystem(V,Vvecs,Vvals)
	invsigma = Vvecs * (Vvecs' :/ Vvals')
	q2 = 0.5 :* rowsum((X*invsigma):*X)
	d1 = -0.5 :* (rows(V)*log(2*pi()) :+ sum(log(Vvals)) )
	return(exp(d1:-q2))
}

// expand a matrix into all permutations
`NM' stmixed_expandmatrix(`NM' x,| `RS' weights )
{
	`RS' nrows, ncols, nexp, index
	`NM' newx
	
	if (rows(x)==1) return(x)
	else {
		nrows = rows(x)
		ncols = cols(x)	
		nexp  = ncols^(nrows-1)
		newx = J(1,nexp,x)
		index = 1
		for (i=1; i<=nrows-1; i++) {
			reps = ncols^(nrows-i)
			xrep = J(reps,1,x[i,])'
			xrep = rowshape(xrep,1)
			newx[i,] = J(1,index,xrep)
			index = index * ncols
		}
		if (args()>1) {
			for (i=2; i<=nrows; i++) newx[1,] = newx[1,] :* newx[i,]
			return(newx[1,])
		}
		else return(newx)
	}	
}

// outer product by column
`NM' stmixed_outerprod_by_col(`NM' x)
{
	`RS' ncols, nrows
	`NM' newx
	
	nrows = rows(x)
	ncols = cols(x)
	newx = J(nrows,ncols^2,.)
	for (i=1;i<=nrows;i++) newx[i,] = rowshape(cross(x[i,],x[i,]),1)
	return(newx)
}

// error check for se(blups) when 0
void stmixed_seblup_check(`RS' x)
{
	if (x==0 | x==.) {
		errprintf("BLUP calculation failed in adaptive quadrature algorithm\n")
		errprintf("Try increasing gh()\n")
		exit(1986)
	}
}

`NM' panelsum12(`NM' x, `NM' info)
{
	`RS' N, cols
	`NM' res
	
	N = rows(info)
	cols = cols(x)
	res = J(N,cols,.)
	for (i=1; i<=N;i++) res[i,] = quadcolsum(x[|info[i,1],.\info[i,2],.|],1)
	return(res)
}

end

exit
