*! version 1.0.0 13may2013 MJC

//AFT models

version 12.1
mata:

real colvector stmixed_aft_na_lnl(pointer(struct stmixed_main scalar) scalar SS)			
{
	struct stmixed_main scalar S
	S = *SS

	Z_dm_nodes = S.random_vars * (cholesky(S.vcv) * S.baseGHnodes)
	//log survival pdf
	if (S.llog) logspdf = (1:/S.g1):*log(S.l_mat :* exp(-Z_dm_nodes)) :+ (1:/S.g1 :- 1):*log(S.t) :- log(S.g1) :- 2:*log(1:+ (S.l_mat:*exp(-Z_dm_nodes):*S.t):^(1:/S.g1))
	else if (S.lnorm) logspdf = -log(S.t:*S.g1:*sqrt(2:*pi())) :- (1:/(2:*(S.g1:^2))):*(log(S.t):-(S.l_mat:+Z_dm_nodes)):^2
	else {
		gg = 1:/(abs(S.kapp):^2)
		z = sign(S.kapp):*(log(S.t):-(S.l_mat:+Z_dm_nodes)):/S.g1
		if (S.kapp==0) logspdf = -log(S.g1:*S.t:*sqrt(2:*pi())) :- (z:^2):/2 
		else logspdf = gg:*log(gg) :- log(S.g1:*S.t:*sqrt(gg):*gamma(gg)) :+ z:*sqrt(gg):- gg:*exp(abs(S.kapp):*z)
	}
	logspdf = S.d :* logspdf

	//log survival
	if (S.llog) logsurv = -log(1 :+ (S.l_mat:*exp(-Z_dm_nodes):*S.t):^(1:/S.g1))
	else if (S.lnorm) logsurv = log(1 :- normal((log(S.t):-(S.l_mat:+Z_dm_nodes)):/S.g1))
	else {
		gg = 1:/(abs(S.kapp):^2)
		z = sign(S.kapp):*(log(S.t):-(S.l_mat:+Z_dm_nodes))
		if (S.kapp>0) logsurv = log(1:-gammap(gg,gg:*exp(abs(S.kapp):*z)))
		else if (S.kapp==0) logsurv = log(1:-normal(z))
		else logsurv = log(gammap(gg,gg:*exp(abs(S.kapp):*z)))
	}
	logsurv = (1:-S.d) :* logsurv

	for (i=1; i<=S.N; i++) SS->jlnodes[i,] = exp(quadcolsum(logspdf[|S.info[i,1],.\S.info[i,2],.|] :+ logsurv[|S.info[i,1],.\S.info[i,2],.|],1))
	return(log((*SS).jlnodes * S.baseGHweights))
}

real colvector stmixed_aft_lnl(pointer(struct stmixed_main scalar) scalar SS)			
{
	struct stmixed_main scalar S
	S = *SS
	
	for (i=1; i<=S.N; i++) {
		
		panelindex 	= (S.info[i,1],.\S.info[i,2],.)
		Z_dm_nodes  = S.random_vars[|panelindex|] * asarray(S.GHNarray,i)

		//log survival pdf		
		if (S.llog) logspdf = (1:/S.g1):*log(S.l_mat[|panelindex|] :* exp(-Z_dm_nodes)) :+ (1:/S.g1 :- 1):*log(S.t[|panelindex|]) :- log(S.g1) :- 2:*log(1:+ (S.l_mat[|panelindex|]:*exp(-Z_dm_nodes):*S.t[|panelindex|]):^(1:/S.g1))
		else if (S.lnorm) logspdf = -log(S.t[|panelindex|]:*S.g1:*sqrt(2:*pi())) :- (1:/(2:*(S.g1:^2))):*(log(S.t[|panelindex|]):-(S.l_mat[|panelindex|]:+Z_dm_nodes)):^2
		else {
			gg = 1:/(abs(S.kapp):^2)
			z = sign(S.kapp):*(log(S.t[|panelindex|]):-(S.l_mat[|panelindex|]:+Z_dm_nodes)):/S.g1
			if (S.kapp==0) logspdf = -log(S.g1:*S.t[|panelindex|]:*sqrt(2:*pi())) :- (z:^2):/2 
			else logspdf = gg:*log(gg) :- log(S.g1:*S.t[|panelindex|]:*sqrt(gg):*gamma(gg)) :+ z:*sqrt(gg):- gg:*exp(abs(S.kapp):*z)
		}
		logspdf = S.d[|panelindex|] :* logspdf

		//log survival
		if (S.llog) logsurv = -log(1 :+ (S.l_mat[|panelindex|]:*exp(-Z_dm_nodes):*S.t[|panelindex|]):^(1:/S.g1))
		else if (S.lnorm) logsurv = log(1 :- normal((log(S.t[|panelindex|]):-(S.l_mat[|panelindex|]:+Z_dm_nodes)):/S.g1))
		else {
			gg = 1:/(abs(S.kapp):^2)
			z = sign(S.kapp):*(log(S.t[|panelindex|]):-(S.l_mat[|panelindex|]:+Z_dm_nodes))
			if (S.kapp>0) logsurv = log(1:-gammap(gg,gg:*exp(abs(S.kapp):*z)))
			else if (S.kapp==0) logsurv = log(1:-normal(z))
			else logsurv = log(gammap(gg,gg:*exp(abs(S.kapp):*z)))
		}
		logsurv = (1:-S.d[|panelindex|]) :* logsurv			
		S.jlnodes[i,] = exp(quadcolsum(logspdf :+ logsurv,1))
	}	
	SS->jlnodes = S.jlnodes :* rowshape(stmixed_dmvnorm(S.stackednodes,S.vcv),S.N) 	
	return(log(S.detcholVCVs :* (*SS).jlnodes  * S.adaptbaseGHweights))
}

end
