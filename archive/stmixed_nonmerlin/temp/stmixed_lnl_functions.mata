*! version 1.0.0 MJC

local TS transmorphic scalar
local RS real scalar
local RC real colvector
local RM real matrix
local PGML pointer(struct GML_struct scalar) scalar
local SGML struct GML_struct scalar

version 12.1
mata:

/*
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
	if ((*S).hasdelentry) {
		`RM' jll
		`RS' lnf2
		jll = (*(*S).refinelnlP)(S)
		S->lnf1 = quadcolsum(jll,1)
		lnf2 = (*(*S).refinelnldelentP)(S)
		S->lnf2 = quadcolsum(lnf2,1)
		jll[(*S).t0idindex,] = jll[(*S).t0idindex,] :- lnf2
		lnf = moptimize_util_sum(M,jll)
	}	
	else S->lnf1 = lnf = moptimize_util_sum(M,(*(*S).refinelnlP)(S))
	if (todo==0) return
}
*/

//================================================================================================================================================================//
// parameter linear predictors

void stmixed_exp_bmat(	`PGML' pGML,
						`TS' M,
						`RR' b)
{
	main_xb = moptimize_util_xb(M,b,1)
	if ((*pGML).Ncoefeq[1]>1) pGML->xb[,1] = main_xb
	else pGML->xb[,1] = J((*pGML).Nobs,1,main_xb)
	
	gml_fillvcv(pGML,M,b,2)
}

void stmixed_weib_bmat(	`PGML' pGML,
						`TS' M,
						`RR' b)
{
	main_xb = moptimize_util_xb(M,b,1)
	if ((*pGML).Ncoefeq[1]>1) pGML->xb[,1] = main_xb
	else pGML->xb[,1] = J((*pGML).Nobs,1,main_xb)
	
	main_xb = exp(moptimize_util_xb(M,b,2))
	if ((*pGML).Ncoefeq[2]>1) pGML->xb[,2] = main_xb
	else pGML->xb[,2] = J((*pGML).Nobs,1,main_xb)

	gml_fillvcv(pGML,M,b,3)
}

void stmixed_gomp_bmat(	`PGML' pGML,
						`TS' M,
						`RR' b)
{
	main_xb = moptimize_util_xb(M,b,1)
	if ((*pGML).Ncoefeq[1]>1) pGML->xb[,1] = main_xb
	else pGML->xb[,1] = J((*pGML).Nobs,1,main_xb)
	
	main_xb = moptimize_util_xb(M,b,2)
	if ((*pGML).Ncoefeq[2]>1) pGML->xb[,2] = main_xb
	else pGML->xb[,2] = J((*pGML).Nobs,1,main_xb)
	
	gml_fillvcv(pGML,M,b,3)
}

void stmixed_fpm_bmat(	`PGML' pGML,
						`TS' M,
						`RR' b)
{
	pGML->xb[,1] = moptimize_util_xb(M,b,1)
	pGML->xb[,2] = moptimize_util_xb(M,b,2)
	pGML->xb[,3] = asarray((*pGML).X,3) * (b[,((*pGML).Ncoefeq[1]+1)..((*pGML).Ncoefeq[1]:+(*pGML).Ncoefeq[2]-1)])'

	//pGML->xb[,2] = moptimize_util_xb(M,b,2)
	/*if ((*pGML).hasdelentry) {
		pGML->xb[,3] = moptimize_util_xb(M,b,3)[(*S).t0index,]
		ind = 4
	}
	else*/
	gml_fillvcv(pGML,M,b,3)
}

void stmixed_llog_lnorm_bmat(	`PGML' pGML,`PS' S, 
								scalar M,
								rowvector b)
{
	main_xb = moptimize_util_xb(M,b,1)
	if ((*pGML).Ncoefeq[1]>1) pGML->xb[,1] = main_xb
	else pGML->xb[,1] = J((*pGML).Nobs,1,main_xb)
	
	main_xb = exp(moptimize_util_xb(M,b,2))
	if ((*pGML).Ncoefeq[2]>1) pGML->xb[,2] = main_xb
	else pGML->xb[,2] = J((*pGML).Nobs,1,main_xb)

	gml_fillvcv(pGML,M,b,3)
}

void stmixed_ggamma_bmat(	`PGML' pGML,
							`TS' M,
							`RR' b)
{
	main_xb = moptimize_util_xb(M,b,1)
	if ((*pGML).Ncoefeq[1]>1) pGML->xb[,1] = main_xb
	else pGML->xb[,1] = J((*pGML).Nobs,1,main_xb)
	
	main_xb = exp(moptimize_util_xb(M,b,2))
	if ((*pGML).Ncoefeq[2]>1) pGML->xb[,2] = main_xb
	else pGML->xb[,2] = J((*pGML).Nobs,1,main_xb)	

	main_xb = moptimize_util_xb(M,b,3)
	if ((*pGML).Ncoefeq[3]>1) pGML->xb[,3] = main_xb
	else pGML->xb[,3] = J((*pGML).Nobs,1,main_xb)	

	gml_fillvcv(pGML,M,b,4)
}


//================================================================================================================================================================//
//log-likelihood at the observations level for all models

`RM' stmixed_weib_ll(`PGML' pGML, `RM' linpred)
{
	`SGML' GML
	GML = (*pGML)
	return(GML.y[,2] :* (linpred :+ log(GML.xb[,2]) :+ (GML.xb[,2]:-1):*log(GML.y[,1]) :+ log(GML.y[,1])) :- exp(linpred):*GML.y[,1]:^(GML.xb[,2]))
}

`RM' stmixed_weib_ll_d_lnlambda(`PGML' pGML, `RM' linpred)
{
	`SGML' GML
	GML = *pGML
	return(GML.y[,2]  :- exp(linpred):*GML.y[,1]:^(GML.xb[,2]))
}

`RM' stmixed_weib_ll_d_lngamma(`PGML' pGML, `RM' linpred)
{
	`SGML' GML
	GML = *pGML
	return(GML.y[,2] :* (1 :+ GML.xb[,2] :* log(GML.y[,1])) :- exp(linpred):*GML.y[,1]:^(GML.xb[,2]) :* log(GML.y[,1]) :* GML.xb[,2])
}

`RM' stmixed_weib_ll_d2_lnlambda(`PGML' pGML, `RM' linpred)
{
	`SGML' GML
	GML = *pGML
	return(-exp(linpred):*GML.y[,1]:^(GML.xb[,2]))
}

`RM' stmixed_weib_ll_d2_lngamma(`PGML' pGML, `RM' linpred)
{
	`SGML' GML
	GML = *pGML
	return(GML.y[,2] :* GML.xb[,2] :* log(GML.y[,1]) :- exp(linpred):*GML.y[,1]:^(GML.xb[,2]) :* log(GML.y[,1]) :* GML.xb[,2] :* (1 :+ GML.xb[,2] :* log(GML.y[,1])))
}

`RM' stmixed_weib_ll_dllamdlgam(`PGML' pGML, `RM' linpred)
{
	`SGML' GML
	GML = *pGML
	return(-exp(linpred):*GML.y[,1]:^(GML.xb[,2]) :* log(GML.y[,1]) :* GML.xb[,2])
}

`RM' stmixed_weib_ll_d_lnsigma(`PGML' pGML, `RM' linpred, `RM' ui)
{
	`SGML' GML
	GML = *pGML
	return(ui :* (GML.y[,2] :- exp(linpred):*GML.y[,1]:^(GML.xb[,2])))
}

`RM' stmixed_weib_ll_d2_lnsigma(`PGML' pGML, `RM' linpred, `RM' ui)
{
	`SGML' GML
	GML = *pGML
	return(ui :* (GML.y[,2] :- exp(linpred):*GML.y[,1]:^(GML.xb[,2])) :- exp(linpred):*GML.y[,1]:^(GML.xb[,2]) :* ui:^2)
}

`RM' stmixed_weib_ll_d2_lnsigma_unstr(`PGML' pGML, `RM' linpred, `RM' ui, `RM' ui2)
{
	`SGML' GML
	GML = *pGML
	return(ui2 :* (GML.y[,2] :- exp(linpred):*GML.y[,1]:^(GML.xb[,2])) :- exp(linpred):*GML.y[,1]:^(GML.xb[,2]) :* ui:^2)
}

`RM' stmixed_weib_ll_dlnl_dlnsig(`PGML' pGML, `RM' linpred, `RM' ui)
{
	`SGML' GML
	GML = *pGML
	return(-exp(linpred):*GML.y[,1]:^(GML.xb[,2]) :* ui)
}

`RM' stmixed_weib_ll_dlnga_dlnsig(`PGML' pGML, `RM' linpred, `RM' ui)
{
	`SGML' GML
	GML = *pGML
	return(-exp(linpred):*GML.y[,1]:^(GML.xb[,2]) :* ui :* log(GML.y[,1]) :* GML.xb[,2])
}

`RM' stmixed_weib_ll_dlnsig_dlnsig(`PGML' pGML, `RM' linpred, `RM' ui, `RM' ui2)
{
	`SGML' GML
	GML = *pGML
	return(-exp(linpred):*GML.y[,1]:^(GML.xb[,2]) :* ui :* ui2)
}

`RM' stmixed_weib_ll_dlns_dlns_unstr(`PGML' pGML, `RM' linpred, `RM' ui, `RM' ui2,| `RM' ui3)
{
	`SGML' GML
	GML = *pGML
	if (args()==4) return((-exp(linpred):*GML.y[,1]:^(GML.xb[,2])) :* ui :* ui2)
	else return((-exp(linpred):*GML.y[,1]:^(GML.xb[,2])) :* ui :* ui2 :+ (GML.y[,2] :- exp(linpred):*GML.y[,1]:^(GML.xb[,2])):* ui3)
}

`RM' stmixed_fpm_ll(`PGML' pGML, `RM' linpred)
{
	`SGML' GML
	GML = *pGML
	/*if (GML.hasbhazard) return(GML.y[,2] :* log(S.bhazard :+ exp(log(GML.xb[,2]):-log(GML.y[,1]) :+ linpred) :- exp(linpred)))
	else*/
	return(GML.y[,2] :* (log(GML.xb[,3]) :+ linpred :+ GML.xb[,2]) :- exp(linpred:+ GML.xb[,2]))
}

`RM' stmixed_fpm_d_xb(`PGML' pGML, `RM' linpred)
{
	`SGML' GML
	GML = *pGML
	return(GML.y[,2] :- exp(linpred :+ GML.xb[,2]))
}

`RM' stmixed_fpm_d_rcs(`PGML' pGML, `RM' linpred, `RS' cind)
{
	`SGML' GML
	GML = *pGML
	if (cind<GML.Ncoefeq[2]) return(GML.y[,2] :* asarray(GML.X,3)[,cind] :/ GML.xb[,3] :+ asarray(GML.X,2)[,cind] :* (GML.y[,2] :- exp(linpred :+ GML.xb[,2])))
	else return(GML.y[,2] :- exp(linpred :+ GML.xb[,2]))
}

`RM' stmixed_fpm_d_lnsigma(`PGML' pGML, `RM' linpred, `RM' ui)
{
	`SGML' GML
	GML = *pGML
	return(ui :* (GML.y[,2] :- exp(linpred :+ GML.xb[,2])))
}

`RM' stmixed_fpm_d2_xb(`PGML' pGML, `RM' linpred)
{
	`SGML' GML
	GML = *pGML
	return(-exp(linpred :+ GML.xb[,2]))
}

`RM' stmixed_fpm_d2_rcs(`PGML' pGML, `RM' linpred, `RS' cind1, `RS' cind2)
{
	`SGML' GML
	GML = *pGML
	if (cind1==GML.Ncoefeq[2] & cind2==GML.Ncoefeq[2]) return(-exp(linpred :+ GML.xb[,2]))
	else if (cind1==GML.Ncoefeq[2] & cind2<GML.Ncoefeq[2]) return(-exp(linpred :+ GML.xb[,2]) :* asarray(GML.X,2)[,cind2])
	else return(-GML.y[,2] :* asarray(GML.X,3)[,cind1]:* asarray(GML.X,3)[,cind2] :/ (GML.xb[,3]:^2) :- asarray(GML.X,2)[,cind1]:* asarray(GML.X,2)[,cind2] :* exp(linpred :+ GML.xb[,2]))
}

`RM' stmixed_fpm_d2_lnsigma(`PGML' pGML, `RM' linpred, `RM' ui)
{
	`SGML' GML
	GML = *pGML
	return(ui :* (GML.y[,2] :- exp(linpred :+ GML.xb[,2])) :- exp(linpred :+ GML.xb[,2]) :* ui:^2)
}

`RM' stmixed_fpm_dxb_drcs(`PGML' pGML, `RM' linpred, `RS' cind1, `RS' cind2)
{
	`SGML' GML
	GML = *pGML
	return(-exp(linpred :+ GML.xb[,2]) :* asarray(GML.X,2)[,cind1] :* asarray(GML.X,1)[,cind2])
}

`RM' stmixed_fpm_d_lnsigma_d_xb(`PGML' pGML, `RM' linpred, `RM' ui)
{
	`SGML' GML
	GML = *pGML
	return(-exp(linpred :+ GML.xb[,2]) :* ui)
}

`RM' stmixed_fpm_d_lnsigma_d_rcs(`PGML' pGML, `RM' linpred, `RM' ui, `RS' cind)
{
	`SGML' GML
	GML = *pGML
	return(-exp(linpred :+ GML.xb[,2]) :* ui :* asarray(GML.X,2)[,cind])
}

`RM' stmixed_fpm_dlnsig_dlnsig(`PGML' pGML, `RM' linpred, `RM' ui, `RM' ui2)
{
	`SGML' GML
	GML = *pGML
	return(-exp(linpred :+ GML.xb[,2]) :* ui :* ui2)
}

/*
numeric matrix stmixed_exp_ll(`PS' pS, `RC' linpred)
{
	`SS' S
	S = *pS
	return(S.d :* (linpred :+ log(S.t)) :- exp(linpred):*S.t)
}

numeric matrix stmixed_gomp_ll(`PS' pS, `RC' linpred)
{
	`SS' S
	S = *pS
	return(S.d :* (linpred :+ S.g1:*S.t :+ log(S.t)) :- exp(linpred):*(1/S.g1):*(exp(S.g1:*S.t):-1))
}

numeric matrix stmixed_llog_ll(`PS' pS, `RC' linpred)
{
	`SS' S
	S = *pS
	return(S.d :* ((1:/S.g1):*(-linpred) :+ (1:/S.g1 :- 1):*log(S.t) :- log(S.g1) :- 2:*log(1:+ (exp(-linpred):*S.t):^(1:/S.g1)))  :- (1:-S.d)  :* log(1 :+ (exp(-linpred):*S.t):^(1:/S.g1)))
}

numeric matrix stmixed_lnorm_ll(`PS' pS, `RC' linpred)
{
	`SS' S
	S = *pS
	return(S.d :* (-log(S.t:*S.g1:*sqrt(2:*pi())) :- (1:/(2:*(S.g1:^2))):*(log(S.t):-linpred):^2) :+ (1:-S.d) :* (log(1 :- normal((log(S.t):-linpred):/S.g1))) )
}

numeric matrix stmixed_ggamma_ll(`PS' pS, `RC' linpred)
{
	`SS' S
	S = *pS
	
	gg = 1:/(abs(S.kapp):^2)
	z = sign(S.kapp):*(log(S.t):-linpred):/S.g1
	if (S.kapp==0) logspdf = -log(S.g1:*S.t:*sqrt(2:*pi())) :- (z:^2):/2 
	else logspdf = gg:*log(gg) :- log(S.g1:*S.t:*sqrt(gg):*gamma(gg)) :+ z:*sqrt(gg):- gg:*exp(abs(S.kapp):*z)
	
	if (S.kapp>0) logsurv = log(1:-gammap(gg,gg:*exp(abs(S.kapp):*z)))
	else if (S.kapp==0) logsurv = log(1:-normal(z))
	else logsurv = log(gammap(gg,gg:*exp(abs(S.kapp):*z)))
	
	return(S.d :* (logspdf) :+ (1:-S.d) :* (logsurv))
}

//================================================================================================================================================================//
//log-survival at the observations level for all models at entry times for delayed entry models

numeric matrix stmixed_weib_lsurv_t0(`PS' pS, `RC' linpred)
{
	`SS' S
	S = *pS
	return(-exp(linpred):*S.t0:^(S.g1))
}

numeric matrix stmixed_exp_lsurv_t0(`PS' pS, `RC' linpred)
{
	`SS' S
	S = *pS
	return(-exp(linpred):*S.t0)
}

numeric matrix stmixed_gomp_lsurv_t0(`PS' pS, `RC' linpred)
{
	`SS' S
	S = *pS
	return(-exp(linpred):*(1/S.g1):*(exp(S.g1:*S.t0):-1))
}

numeric matrix stmixed_fpm_lsurv_t0(`PS' pS, `RC' linpred)
{
	`SS' S
	S = *pS
	return(-exp(linpred))
}

numeric matrix stmixed_llog_lsurv_t0(`PS' pS, `RC' linpred)
{
	`SS' S
	S = *pS
	return(log(1 :+ (exp(-linpred):*S.t0):^(1:/S.g1)))
}

numeric matrix stmixed_lnorm_lsurv_t0(`PS' pS, `RC' linpred)
{
	`SS' S
	S = *pS
	return(-log(1 :- normal((log(S.t0):-linpred):/S.g1)))
}

numeric matrix stmixed_ggamma_lsurv_t0(`PS' pS, `RC' linpred)
{
	`SS' S
	S = *pS
	
	gg = 1:/(abs(S.kapp):^2)
	z = sign(S.kapp):*(log(S.t0):-linpred):/S.g1
	if (S.kapp>0) return(-log(1:-gammap(gg,gg:*exp(abs(S.kapp):*z))))
	else if (S.kapp==0) return(-log(1:-normal(z)))
	else return(-log(gammap(gg,gg:*exp(abs(S.kapp):*z))))
}
*/

end
