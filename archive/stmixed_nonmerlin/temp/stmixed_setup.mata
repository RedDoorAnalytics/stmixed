*! version 1.0.0 13may2013 MJC

version 12.1

local SS 	string scalar
local RS	real scalar
local RM	real matrix
local RC	real colvector
local TR	transmorphic
local ps	pointer scalar
local RR 	real rowvector
local PGML	pointer(struct GML_struct scalar) scalar

mata:

void stmixed_setup()
{
	`SS' touse
	`PGML' pGML
	
	stata("tempname gml_struct")
	rmexternal(st_local("gml_struct"))
	pGML = crexternal(st_local("gml_struct"))
	touse = st_local("touse")
	smodel = st_local("smodel")
	
	GML = gml_init()
	gml_Nobs(GML,st_numscalar(st_local("obs")))
	gml_Nlevels(GML,strtoreal(st_local("Nlevels")))
	gml_panelinfo(GML,"levelvar","touse")	
	gml_Nres_at_levels(GML,strtoreal(tokens(st_local("Nres")))')
	gml_init_gauss_hermite(GML,strtoreal(st_local("gh")))
	gml_P_lnl(GML,stmixed_lnl_pointer(smodel))
	gml_init_lnli(GML)
	gml_P_betas(GML,stmixed_fillbeta_pointer(smodel))
	gml_vcv_mats(GML)
	gml_vcv_structures(GML,st_local("covs"))
	gml_qind(GML)
	gml_init_Z(GML,"varnames_",touse)
	gml_adaptiveGH(GML,st_local("nonadapt")=="")
	gml_todo(GML,strtoreal(st_local("todo")))
	gml_NmleqnsFEs(GML,strtoreal(st_local("NmleqnsFE")))	//!!generalise
	gml_init_X(GML,"Xdm_",touse)
	gml_P_scoreFE(GML,stmixed_get_PscoreFE(smodel))
	gml_P_scoreSigma(GML,stmixed_get_PscoreSigma(smodel))
	gml_chainrule(GML)
	gml_init_G(GML)
	if (strtoreal(st_local("todo"))==2) gml_P_hessian(GML,hessianPmat(smodel,strtoreal(tokens(st_local("Nres")))'))
	gml_bhazard(GML,st_local("bhazard")!="")
	
	gml_init_xb(GML)
	gml_init_y(GML,(st_data(.,"_t",touse),st_data(.,"_d",touse)))
	
	if (smodel=="fpm") {
		gml_chainrule(GML,(1\0))
		gml_NmleqnsFEs(GML,strtoreal(st_local("NmleqnsFE"))-1)
	}
	swap(*pGML,GML)
	
		
	//Survival model info 
	/*S.smodel = st_local("smodel")
	S.hasbhazard = st_local("bhazard")!=""
	if (S.hasbhazard) st_view(S.bhazard,.,st_local("bhazard"),touse)
	S.hascovs = st_local("varnames_0")!=""		//need for ml equation (scalar or vector)
	S.hasdelentry = strtoreal(st_local("delentry"))
	if (S.hasdelentry) {
		S.t0 = st_data(.,"_t0",st_local("t0touse"))
		S.t0index = st_data(.,st_local("t0index"),st_local("t0touse"))
		S.t0idindex = uniqrows(st_data(.,st_local("levelvar1"),st_local("t0touse")))
		S.Nt0obs = strtoreal(st_local("Nt0obs"))
		S.fpm = S.smodel=="fpm"
	}*/
	
	//if (S.hasdelentry) S.lnf2 = J(1,1,.)
	
	//Initialise FE parameter vectors
	/*if (S.smodel=="fpm") {
		S.dxb = J((*pGML).Nobs,1,.)
		if (S.hasdelentry) S.xb0 = J((*pGML).Nobs,1,.)
	}
	if (S.smodel!="e" & S.smodel!="fpm") S.g1 = J(1,1,.)
	if (S.smodel=="gamma") S.kapp = J(1,1,.)
	*/
	//GH setup
	//stmixed_gh_setup(&S,Ngh,S.adapt,touse)
	/*if (S.adapt) {
		S.iter = 0
		S.adaptit = strtoreal(st_local("adaptit"))
		S.atol = strtoreal(st_local("retolerance"))
		S.showadapt = st_local("showadapt")!=""
	}*/
	
	//if (S.hasdelentry) S.plogs 		= stmixed_lnl_delent_pointer(S.smodel)
	
	//if (S.getblups) stmixed_getblups(p)
}

/*
void stmixed_gh_setup(	`PS' S,
						real scalar NGHnodes,
						real scalar quadtype,
						string scalar touse)
{
	`RM' qmat, x, w, baseGHnodes, baseGHweights
	`TR' aghnodes
	`RC' Nres
	
	Nres = (*S).GML.Nres
	S->baseGHnodes = S->baseGHweights = asarray_create("real",1)
	
	qmat = _gauss_hermite_nodes(NGHnodes)
	for (i=1; i<=(*S).GML.Nlevels; i++) {
		x = J(Nres[i],1,qmat[1,]):*sqrt(2)
		baseGHnodes = expand_matrix(x)
		asarray(S->baseGHnodes,i,baseGHnodes)
		
		w = J(Nres[i],1,qmat[2,]):/sqrt(pi())
		baseGHweights = expand_matrix(w,1)'
		asarray(S->baseGHweights,i,baseGHweights)
	}
	
	//Adaptive section
	/*if (quadtype) {
		aghnodes 	= asarray_create("real",1)
		S->detcholVCVs = J((*S).N,1,1)
		for (i=1; i<=(*S).N; i++) asarray(aghnodes,i,baseGHnodes)
		S->GHNarray = aghnodes
		S->adaptbaseGHweights = ((2:*pi()):^((*S).Nres:/2):*exp(quadcolsum(baseGHnodes:^2):/2) :* (*S).baseGHweights')'
		S->stackednodes = J((*S).N,1,(*S).baseGHnodes')
		
		if ((*S).hasdelentry) {
			`TR' del_aghnodes
			del_aghnodes 	= asarray_create("real",1)
			for (i=1; i<=(*S).Ndelentid; i++) asarray(del_aghnodes,i,baseGHnodes)
			S->del_GHNarray = del_aghnodes
			S->del_detcholVCVs = J((*S).Ndelentid,1,1)
			S->del_stackednodes = J((*S).Ndelentid,1,baseGHnodes')
			S->del_adaptbaseGHweights = (*S).adaptbaseGHweights
		}
	}*/
}
*/
pointer scalar stmixed_fillbeta_pointer(`SS' smodel) 
{
	if (smodel=="w") 	return(&stmixed_weib_bmat())	
	if (smodel=="gom") 	return(&stmixed_gomp_bmat())	
	if (smodel=="e") 	return(&stmixed_exp_bmat())
	if (smodel=="fpm")	return(&stmixed_fpm_bmat())
	if (smodel=="llogistic") return(&stmixed_llog_lnorm_bmat())
	if (smodel=="lnormal") return(&stmixed_llog_lnorm_bmat())
	if (smodel=="gamma") return(&stmixed_ggamma_bmat())
}

pointer scalar stmixed_lnl_pointer(`SS' smodel)
{
	if (smodel=="fpm") return(&stmixed_fpm_ll())	
	else if (smodel=="llogistic") return(&stmixed_llog_ll())
	else if (smodel=="lnormal") return(&stmixed_lnorm_ll())
	else if (smodel=="gamma") return(&stmixed_ggamma_ll())
	else if (smodel=="w") return(&stmixed_weib_ll())
	else if (smodel=="e") return(&stmixed_exp_ll())
	else return(&stmixed_gomp_ll())
}

pointer colvector stmixed_get_PscoreFE(`SS' smodel)
{
	if (smodel=="fpm") return(&stmixed_fpm_d_xb()\&stmixed_fpm_d_rcs())	
	else if (smodel=="llogistic") return(&stmixed_llog_ll())
	else if (smodel=="lnormal") return(&stmixed_lnorm_ll())
	else if (smodel=="gamma") return(&stmixed_ggamma_ll())
	else if (smodel=="w") {
		return(&stmixed_weib_ll_d_lnlambda()\&stmixed_weib_ll_d_lngamma())
	}
	else if (smodel=="e") return(&stmixed_exp_ll())
	else return(&stmixed_gomp_ll())
}

pointer colvector stmixed_get_PscoreSigma(`SS' smodel)
{
	if (smodel=="fpm") return(&stmixed_fpm_d_lnsigma())	
	else if (smodel=="llogistic") return(&stmixed_llog_ll())
	else if (smodel=="lnormal") return(&stmixed_lnorm_ll())
	else if (smodel=="gamma") return(&stmixed_ggamma_ll())
	else if (smodel=="w") {
		return(&stmixed_weib_ll_d_lnsigma())
	}
	else if (smodel=="e") return(&stmixed_exp_ll())
	else return(&stmixed_gomp_ll())
}

pointer scalar stmixed_lnl_delent_pointer(`SS' smodel)
{
	if (smodel=="fpm") return(&stmixed_fpm_lsurv_t0())	
	else if (smodel=="llogistic") return(&stmixed_llog_lsurv_t0())
	else if (smodel=="lnormal") return(&stmixed_lnorm_lsurv_t0())
	else if (smodel=="gamma") return(&stmixed_ggamma_lsurv_t0())
	else if (smodel=="w") return(&stmixed_weib_lsurv_t0())
	else if (smodel=="e") return(&stmixed_exp_lsurv_t0())
	else return(&stmixed_gomp_lsurv_t0())
}

`TR' hessianPmat(`SS' smodel, `RC' Nres) 
{
	Nlevels = rows(Nres)
	covs = st_local("covs")
	covariance = tokens(covs):=="independent"\tokens(covs):=="exchangeable"\tokens(covs):=="unstructured"

	if (smodel=="w") {
		hessian = &stmixed_weib_ll_d2_lnlambda(),&stmixed_weib_ll_dllamdlgam()\&stmixed_weib_ll_dllamdlgam(),&stmixed_weib_ll_d2_lngamma()
		for (i=1;i<=Nlevels;i++) {
			if (covariance[1,i]) {
				for (j=1;j<=Nres[i];j++) {
					cov = &stmixed_weib_ll_dlnl_dlnsig(),&stmixed_weib_ll_dlnga_dlnsig()
					
					//higher levels
					lev = 1
					while (lev<i) {
						if (covariance[1,lev]) {
							for (re=1;re<=Nres[lev];re++) {
								cov = cov,&stmixed_weib_ll_dlnsig_dlnsig()
							}
						}
						else if (covariance[2,lev]) {
						
						}
						else if (covariance[3,lev]) {
						
						}
						else cov = cov,&stmixed_weib_ll_dlnsig_dlnsig()
						lev++
					}
					reind = 1
					while (reind<j) {
						cov = cov,&stmixed_weib_ll_dlnsig_dlnsig()
						reind++
					}
					hessian = hessian,cov'\cov,&stmixed_weib_ll_d2_lnsigma()
				}
			
			}
			else if (covariance[2,i]) {
			
			}
			else if (covariance[3,i]) {
			
				Hresind = 1
				for (j=1;j<=Nres[i];j++) {
					
					k=1
					while (k<=j) {
					
						cov = &stmixed_weib_ll_dlnl_dlnsig(),&stmixed_weib_ll_dlnga_dlnsig()
						
						//higher levels
						lev = 1
						while (lev<i) {
							if (covariance[1,lev]) {
								for (re=1;re<=Nres[lev];re++) {
									cov = cov,&stmixed_weib_ll_dlnsig_dlnsig()
								}
							}
							else if (covariance[2,lev]) {
							
							}
							else if (covariance[3,lev]) {
								for (re=1;re<=Nres[lev];re++) {
									reind = 1
									while (reind<j) {
										for (reind2=1;reind2<=reind;reind2++) cov = cov,&stmixed_weib_ll_dlnsig_dlnsig()
										reind++
									}
								}
							}
							else cov = cov,&stmixed_weib_ll_dlnsig_dlnsig()
							lev++
						}
						//same level
						for (reind2=1;reind2<Hresind;reind2++) cov = cov,&stmixed_weib_ll_dlns_dlns_unstr()
						/*reind = 1
						while (reind<j) {
							for (reind2=1;reind2<=reind;reind2++) cov = cov,&stmixed_weib_ll_d_lnsigma()
							reind++
						}
						*/
						hessian = hessian,cov'\cov,&stmixed_weib_ll_d2_lnsigma_unstr()
						k++
						
						Hresind++
					}
				}
			
			
			}
			else {
				cov = &stmixed_weib_ll_dlnl_dlnsig(),&stmixed_weib_ll_dlnga_dlnsig()
				//higher levels
				lev = 1
				while (lev<i) {
					if (covariance[1,lev]) {
						for (re=1;re<=Nres[lev];re++) {
							cov = cov,&stmixed_weib_ll_dlnsig_dlnsig()
						}
					}
					else if (covariance[2,lev]) {
					
					}
					else if (covariance[3,lev]) {
					
					}
					else cov = cov,&stmixed_weib_ll_dlnsig_dlnsig()
					lev++
				}
				hessian = hessian,cov'\cov,&stmixed_weib_ll_d2_lnsigma()
			
			}
		}
	}
	if (smodel=="fpm") {
		hessian = &stmixed_fpm_d2_xb(),&stmixed_fpm_dxb_drcs()\&stmixed_fpm_dxb_drcs(),&stmixed_fpm_d2_rcs()
		for (i=1;i<=Nlevels;i++) {
			if (covariance[1,i]) {
				for (j=1;j<=Nres[i];j++) {
					cov = &stmixed_fpm_d_lnsigma_d_xb(),&stmixed_fpm_d_lnsigma_d_rcs()
					
					//higher levels
					lev = 1
					while (lev<i) {
						if (covariance[1,lev]) {
							for (re=1;re<=Nres[lev];re++) {
								cov = cov,&stmixed_fpm_dlnsig_dlnsig()
							}
						}
						else if (covariance[2,lev]) {
						
						}
						else if (covariance[3,lev]) {
						
						}
						else cov = cov,&stmixed_fpm_dlnsig_dlnsig()
						lev++
					}
					reind = 1
					while (reind<j) {
						cov = cov,&stmixed_fpm_dlnsig_dlnsig()
						reind++
					}
					hessian = hessian,cov'\cov,&stmixed_fpm_d2_lnsigma()
				}
			
			}
			else if (covariance[2,i]) {
			
			}
			else if (covariance[3,i]) {
			
			}
			else {
				cov = &stmixed_fpm_d_lnsigma_d_xb(),&stmixed_fpm_d_lnsigma_d_rcs()
				//higher levels
				lev = 1
				while (lev<i) {
					if (covariance[1,lev]) {
						for (re=1;re<=Nres[lev];re++) {
							cov = cov,&stmixed_fpm_dlnsig_dlnsig()
						}
					}
					else if (covariance[2,lev]) {
					
					}
					else if (covariance[3,lev]) {
					
					}
					else cov = cov,&&stmixed_fpm_dlnsig_dlnsig()
					lev++
				}
				hessian = hessian,cov'\cov,&stmixed_fpm_d2_lnsigma()
			
			}
		}
	}
	return(hessian)
}

function myexit()
{
	printf("\nExited normally\n\n")
	exit(1986)
}

end

exit

