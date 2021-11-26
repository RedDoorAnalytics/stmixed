*! version 1.0.0 09may2013 MJC

//FPM likelihood evaluator functions

version 12.1
mata:

real colvector stmixed_fpm_na_lnl(pointer(struct stmixed_main scalar) scalar SS)			
{
	struct stmixed_main scalar S
	S = *SS

	Z_dm_nodes = S.random_vars * (cholesky(S.vcv) * S.baseGHnodes)
	if (S.hasbhazard) loghaz = S.d :* log(S.bhazard :+ exp(log(S.dxb):-log(S.t) :+ S.l_mat :+ Z_dm_nodes))
	else loghaz = S.d :* (log(S.dxb) :+ S.l_mat :+ Z_dm_nodes)
	logsurv = -exp(S.l_mat :+ Z_dm_nodes)
	SS->jlnodes[,] = exp(panelsum(loghaz :+ logsurv,S.info))
	return(log((*SS).jlnodes * S.baseGHweights))
}

real colvector stmixed_fpm_lnl(pointer(struct stmixed_main scalar) scalar SS)			
{
	struct stmixed_main scalar S
	S = *SS
	
	for(i=1;i<=S.N;i++) {
		panelindex = (S.info[i,1],.\S.info[i,2],.)
		Z_dm_nodes = S.random_vars[|panelindex|] * asarray(S.GHNarray,i)
		if (S.hasbhazard) loghaz = S.d[|panelindex|] :* log(S.bhazard[|panelindex|] :+ exp(log(S.dxb[|panelindex|]) :- log(S.t[|panelindex|]) :+ S.l_mat[|panelindex|] :+ Z_dm_nodes))
		else loghaz = S.d[|panelindex|] :* (log(S.dxb[|panelindex|]) :+ S.l_mat[|panelindex|] :+ Z_dm_nodes)
		logsurv = -exp(S.l_mat[|panelindex|] :+ Z_dm_nodes)
		S.jlnodes[i,] = exp(quadcolsum(loghaz :+ logsurv,1))
	}
	SS->jlnodes = S.jlnodes :* rowshape(stmixed_dmvnorm(S.stackednodes,S.vcv),S.N) 	
	return(log(S.detcholVCVs :* (*SS).jlnodes  * S.adaptbaseGHweights))
}

end
