*! version 1.0.0 ?????2013 MJC

//stmixed prediction Mata commands

/*
History
MJC ?????2013 version 1.0.0
*/

version 12.1

local NM	numeric matrix
local RS	real scalar

mata:

//fill up struct with fitted values then get blups
void stmixed_getblups(pointer(struct stmixed_main scalar) scalar S)
{
	`RS' a, lnf1, lnf2, postseblups
	
	postseblups = st_local("reise")!=""
	final_lnl = st_numscalar("e(ll)")

	S->vcv = st_matrix("e(vcv)")
		
	//Fill up parameter matrices etc.
	if ((*S).smodel=="w" | (*S).smodel=="gom" | (*S).smodel=="e"){
		stata("tempvar gb_lmat")
		stata("qui _predict double "+st_local("gb_lmat")+" if "+st_local("touse")+",xb eq(ln_lambda)")
		S->l_mat = exp(st_data(.,st_local("gb_lmat"),st_local("survtouse")))
		if ((*S).smodel!="e") {
			if ((*S).smodel=="w") {
				stata("local gb_g = [ln_p][_cons]")
				S->g1 = exp(strtoreal(st_local("gb_g")))
			}
			else if ((*S).smodel=="gom") {
				stata("local gb_g = [gamma][_cons]")
				S->g1 = strtoreal(st_local("gb_g"))
			}
		}
		
	}
	if ((*S).smodel=="fpm") {
		stata("tempvar gb_xb")
		stata("qui _predict double "+st_local("gb_xb")+" if "+st_local("touse")+",xb eq(xb)")
		S->l_mat 	= st_data(.,st_local("gb_xb"),st_local("survtouse"))
		stata("tempvar gb_dxb")
		stata("qui _predict double "+st_local("gb_dxb")+" if "+st_local("touse")+",xb eq(dxb)")
		S->dxb 	= st_data(.,st_local("gb_dxb"),st_local("survtouse"))
	}
	
	//adaptive or non-adaptive
	if ((*S).adapt) {
		lnf1 = quadcolsum((*(*S).lnlP)(S),1)
		stmixed_update_quad(S)
		lnf2 = quadcolsum((*(*S).lnlP)(S),1)
		while(abs(lnf1:-lnf2)>(*S).atol) {
			swap(lnf1,lnf2)
			stmixed_update_quad(S)
			lnf2 = quadcolsum((*(*S).lnlP)(S),1)
		}
	}
	else lnf2 = quadcolsum((*(*S).lnlP)(S),1)
	
	if (mreldif(final_lnl,lnf2)>1E-08) {
		errprintf("BLUP calculation failed; estimation data have changed\n")
		exit(459)
	}
	stmixed_post_blups(S,postseblups)	
}

// post blups to Stata
void stmixed_post_blups(pointer(struct stmixed_main scalar) scalar SS, `RS' postseblups)
{
	`RS' Nres, N, csurv, csurvid
	`NM' nodes_i, newweights_comp1, vcv_comp1, test1, vcv_new, blups
	
	N = (*SS).N
	Nres = (*SS).Nres
	blups = J(N,Nres,.)
		
	if ((*SS).adapt) {
		like_j = (*SS).detcholVCVs :* (*SS).jlnodes * (*SS).baseGHweights
		for(i=1;i<=N;i++) {
			nodes_i = asarray((*SS).GHNarray,i)		
			test1 = (*SS).detcholVCVs[i,] :* (*SS).jlnodes[i,]:/like_j[i,]
			blups[i,] = test1 * (nodes_i :* (*SS).baseGHweights')'
			if (postseblups) {
				vcv_new = (test1 * (stmixed_outerprod_by_col(nodes_i') :* (*SS).baseGHweights)) :- rowshape(blups[i,]'*blups[i,],1)
				blups[i,] = diagonal(cholesky(rowshape(vcv_new,Nres)))'
			}
		}	
	}
	else {
		like_j = (*SS).jlnodes * (*SS).baseGHweights
		nodes = cholesky((*SS).vcv)*(*SS).baseGHnodes
		for(i=1;i<=N;i++) {
			test1 = (*SS).jlnodes[i,]:/like_j[i,]
			blups[i,] = test1 * (nodes :* (*SS).baseGHweights')'
			if (postseblups) {
				vcv_new = (test1 * (stmixed_outerprod_by_col(nodes') :* (*SS).baseGHweights)) :- rowshape(blups[i,]'*blups[i,],1)
				stmixed_seblup_check(min(diagonal(vcv_new)))	//error check
				blups[i,] = diagonal(cholesky(rowshape(vcv_new,Nres)))'
			}
		}
	}
	st_store(.,tokens(st_local("getblups")),st_local("posttouse"),blups)
}


end

exit
