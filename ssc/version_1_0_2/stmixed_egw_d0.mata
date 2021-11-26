*! version 1.0.0 MJC

version 12.1
mata:

real colvector stmixed_egw_na_lnl(pointer(struct stmixed_main scalar) scalar SS)			
{
	struct stmixed_main scalar S
	S = *SS

	Z_dm_nodes = S.random_vars * (cholesky(S.vcv) * S.baseGHnodes)
	if (S.notexp) {
		if (S.weib) {
			if (S.hasbhazard) loghaz = S.d :* log(S.bhazard :+ exp(log(S.l_mat) :+ log(S.g1) :+ (S.g1-1):*log(S.t) :+ Z_dm_nodes :+ log(S.t)))
			else loghaz = S.d :* (log(S.l_mat) :+ log(S.g1) :+ (S.g1-1):*log(S.t) :+ Z_dm_nodes :+ log(S.t))
			logsurv = -S.l_mat:*S.t:^(S.g1):*exp(Z_dm_nodes)
		}
		else {
			if (S.hasbhazard) loghaz = S.d :* log(S.bhazard :+ exp(log(S.l_mat) :+ S.g1:*S.t :+ Z_dm_nodes :+ log(S.t)))
			else loghaz = S.d :* (log(S.l_mat) :+ S.g1:*S.t :+ Z_dm_nodes :+ log(S.t))
			logsurv = -S.l_mat:*exp(Z_dm_nodes):*(1/S.g1):*(exp(S.g1:*S.t):-1)
		}
	}
	else {
		if (S.hasbhazard) loghaz = S.d :* log(S.bhazard :+ exp(log(S.l_mat) :+ Z_dm_nodes :+ log(S.t)))
		else loghaz = S.d :* (log(S.l_mat) :+ Z_dm_nodes :+ log(S.t))
		logsurv = -S.l_mat:*S.t:*exp(Z_dm_nodes)
	}
	SS->jlnodes[,] = exp(panelsum(loghaz :+ logsurv,S.info))	
	return(log((*SS).jlnodes * S.baseGHweights))
}

real colvector stmixed_egw_lnl(pointer(struct stmixed_main scalar) scalar SS)			
{
	struct stmixed_main scalar S
	S = *SS
	
	for(i=1;i<=S.N;i++) {
		panelindex = (S.info[i,1],.\S.info[i,2],.)
		Z_dm_nodes = S.random_vars[|panelindex|] * asarray(S.GHNarray,i)
		if (S.notexp) {
			if (S.weib) {
				if (S.hasbhazard) loghaz = S.d[|panelindex|] :* log(S.bhazard[|panelindex|] :+ exp(log(S.l_mat[|panelindex|]) :+ log(S.g1) :+ (S.g1-1):*log(S.t[|panelindex|]) :+ Z_dm_nodes :+ log(S.t[|panelindex|])))
				else loghaz = S.d[|panelindex|] :* (log(S.l_mat[|panelindex|]) :+ log(S.g1) :+ (S.g1-1):*log(S.t[|panelindex|]) :+ Z_dm_nodes :+ log(S.t[|panelindex|]))
				logsurv = -S.l_mat[|panelindex|]:*S.t[|panelindex|]:^(S.g1):*exp(Z_dm_nodes)
			}
			else {
				if (S.hasbhazard) loghaz = S.d[|panelindex|] :* log(S.bhazard[|panelindex|] :+ exp(log(S.l_mat[|panelindex|]) :+ S.g1:*S.t[|panelindex|] :+ Z_dm_nodes :+ log(S.t[|panelindex|])))
				else loghaz = S.d[|panelindex|] :* (log(S.l_mat[|panelindex|]) :+ S.g1:*S.t[|panelindex|] :+ Z_dm_nodes :+ log(S.t[|panelindex|]))
				logsurv = -S.l_mat[|panelindex|]:*exp(Z_dm_nodes):*(1/S.g1):*(exp(S.g1:*S.t[|panelindex|]):-1)
			}
		}
		else {
			if (S.hasbhazard) loghaz = S.d[|panelindex|] :* log(S.bhazard[|panelindex|] :+ exp(log(S.l_mat[|panelindex|]) :+ Z_dm_nodes :+ log(S.t[|panelindex|])))
			else loghaz = S.d[|panelindex|] :* (log(S.l_mat[|panelindex|]) :+ Z_dm_nodes :+ log(S.t[|panelindex|]))
			logsurv = -S.l_mat[|panelindex|]:*S.t[|panelindex|]:*exp(Z_dm_nodes)
		}
		S.jlnodes[i,] = exp(quadcolsum(loghaz :+ logsurv,1))
	}
	SS->jlnodes = S.jlnodes :* rowshape(stmixed_dmvnorm(S.stackednodes,S.vcv),S.N) 	
	return(log(S.detcholVCVs :* (*SS).jlnodes  * S.adaptbaseGHweights))
}

end
