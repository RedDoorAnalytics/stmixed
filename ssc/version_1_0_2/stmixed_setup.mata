*! version 1.0.0 13may2013 MJC

version 12.1

local SS 	string scalar
local RS	real scalar
local NM	numeric matrix

mata:

struct stmixed_main {
	real colvector t, t0, d, l_mat, dxb, xb0, bhazard
	real rowvector xb_fixed
	`NM' random_vars, baseGHnodes, baseGHweights, adaptbaseGHweights, vcv, info, jlnodes, detcholVCVs, stackednodes
	`RS' Nres, ndim, g1, Nobs, N, iter, adaptit, atol, showadapt, covind, covunstr, covexch, hascovs, notexp, weib
	`SS' smodel
	transmorphic GHNarray
	pointer scalar lnlP, fillbetaP, refinelnlP
	`RS' llog, gam, lnorm, kapp
	`RS' getblups, adapt, hasbhazard, refine, reiterate
}

void stmixed_setup()
{
	struct stmixed_main scalar S
	pointer scalar p

	//get tempname and create struct
	stata("tempname stmixed_struct")
	rmexternal(st_local("stmixed_struct"))
	p = crexternal(st_local("stmixed_struct"))
	S.getblups = getblups = st_local("getblups")!=""

	refine = strtoreal(st_local("refine"))
	S.reiterate = strtoreal(st_local("reiterate"))
	touse = st_local("touse")	
	ngh = strtoreal(st_local("gh"))
	S.adapt = st_local("nonadapt")==""
	S.Nres = strtoreal(st_local("nres"))
	S.ndim = ngh:^S.Nres
	S.smodel = st_local("smodel")
	S.Nobs = strtoreal(st_local("nobs"))
	S.N = strtoreal(st_local("N"))
	S.hascovs = st_local("fixed_vars")!=""
	S.hasbhazard = st_local("bhazard")!=""
	aft = strtoreal(st_local("aft"))
	if (aft) {
		S.llog = S.smodel=="llogistic"
		S.lnorm = S.smodel=="lnormal"
		S.gam = S.smodel=="gamma"
	}
	else {
		S.notexp = S.smodel!="e"
		S.weib = S.smodel=="w"
	}
	S.covind 	= st_local("cov")=="ind"
	S.covunstr 	= st_local("cov")=="unstr"
	S.covexch 	= st_local("cov")=="exch"
	
	S.hasbhazard = st_local("bhazard")!=""
	if (S.hasbhazard) st_view(S.bhazard,.,st_local("bhazard"),touse)
	
	st_view(S.t,.,"_t",touse)
	st_view(S.d,.,"_d",touse)
		
	st_view(S.random_vars,.,tokens(st_local("random_vars")),touse)
		
	S.info = panelsetup(st_data(.,st_local("id"),touse), 1)	

	//Initialise parameter vectors
	S.l_mat = J(S.Nobs,1,.)
	if (S.smodel=="fpm") S.dxb = J(S.Nobs,1,.)
	if (S.weib | S.smodel=="gom" | aft) S.g1 = J(1,1,.)
	if (S.gam) S.kapp = J(1,1,.)
	S.vcv = J(S.Nres,S.Nres,.)
	S.jlnodes = J(S.N,S.ndim,.)

	//GH setup
	stmixed_gh_setup(&S,ngh,S.adapt,touse)
	S.iter = 0
	if (S.adapt) S.adaptit = strtoreal(st_local("adaptit"))
	S.atol = strtoreal(st_local("retolerance"))
	S.showadapt = st_local("showadapt")!=""
	
	//pointer to lnl evaluator function
	S.fillbetaP					= stmixed_fillbeta_pointer(S.smodel)
	S.lnlP	 					= stmixed_lnl_pointer(S.smodel,S.adapt)
	if (refine) S.refinelnlP 	= stmixed_lnl_pointer(S.smodel,0)

	//Done 	
	swap((*p), S)
	if (getblups) stmixed_getblups(p)
}

void stmixed_gh_setup(	pointer(struct stmixed_main scalar) scalar S,
						real scalar NGHnodes,
						real scalar quadtype,
						string scalar touse)
{
	qmat = _gauss_hermite_nodes(NGHnodes)
	x = J((*S).Nres,1,qmat[1,]):*sqrt(2)
	S->baseGHnodes   = stmixed_expandmatrix(x)			
	w = J((*S).Nres,1,qmat[2,]):/sqrt(pi())
	S->baseGHweights = stmixed_expandmatrix(w,1)'
	S->ndim = NGHnodes:^(*S).Nres												// No. of permutations for node sequences
	
	//Adaptive section
	if (quadtype) {
		newmu_i = J((*S).N,(*S).Nres,0)
		newtau_i = J((*S).N,(*S).Nres,1)
		aghnodes 	= asarray_create("real",1)
		decomp_i 	= J((*S).Nres,(*S).Nres,0)
		S->detcholVCVs = J((*S).N,1,.)
		
		for (i=1; i<=(*S).N; i++) {
			shift = newmu_i[i,]'
			scale = newtau_i[i,]	
			for(j=1;j<=(*S).Nres;j++) decomp_i[j,j] = scale[1,j]
			nodes_i = shift :+ decomp_i * (*S).baseGHnodes
			S->detcholVCVs[i,] = sqrt(det(decomp_i*decomp_i))
			asarray(aghnodes,i,nodes_i)
			//stack node matrices
			if (i==1) testmat = nodes_i'
			else testmat = testmat\nodes_i'
		}
		S->GHNarray = aghnodes
		S->adaptbaseGHweights = ((2:*pi()):^((*S).Nres:/2):*exp(quadcolsum((*S).baseGHnodes:^2):/2) :* (*S).baseGHweights')'
		S->stackednodes = testmat
	}	
}

pointer scalar stmixed_lnl_pointer(string scalar smodel,real scalar quadtype)
{
 	if (quadtype) {
		if (smodel=="fpm") return(&stmixed_fpm_lnl())
		else if (smodel=="llogistic" | smodel=="lnormal" | smodel=="gamma") return(&stmixed_aft_lnl())
		else return(&stmixed_egw_lnl())	
	}
	else {
		if (smodel=="fpm") return(&stmixed_fpm_na_lnl())	
		else if (smodel=="llogistic" | smodel=="lnormal" | smodel=="gamma") return(&stmixed_aft_na_lnl())
		else return(&stmixed_egw_na_lnl())
	}
}

pointer scalar stmixed_fillbeta_pointer(string scalar smodel)
{
	if (smodel=="w") 	return(&stmixed_w_bmat())	
	if (smodel=="gom") 	return(&stmixed_g_bmat())	
	if (smodel=="e") 	return(&stmixed_e_bmat())
	if (smodel=="fpm")	return(&stmixed_fpm_bmat())
	if (smodel=="llogistic" | smodel=="lnormal" | smodel=="gamma") return(&stmixed_aft_bmat())
}

end

exit
