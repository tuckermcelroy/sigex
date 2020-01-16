sigex.constrainreg <- function(mdl,data.ts,regindex,combos)
{
  
  ##########################################################################
  #
  #	sigex.constrainreg
  # 	    Copyright (C) 2020  Tucker McElroy
  #
  #    This program is free software: you can redistribute it and/or modify
  #    it under the terms of the GNU General Public License as published by
  #    the Free Software Foundation, either version 3 of the License, or
  #    (at your option) any later version.
  #
  #    This program is distributed in the hope that it will be useful,
  #    but WITHOUT ANY WARRANTY; without even the implied warranty of
  #    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  #    GNU General Public License for more details.
  #
  #    You should have received a copy of the GNU General Public License
  #    along with this program.  If not, see <https://www.gnu.org/licenses/>.
  #
  ############################################################################
  
  ################# Documentation #####################################
  #
  #	Purpose: determines constraint matrix for regressors
  #	Inputs:
  #		mdl: the specified sigex model, a list object
  #		data.ts: a T x N matrix ts object
  #   regindex: a list object of N elements, containing indices J_i 
  #     for i=1,...,N of the regressors for the ith series involved.
  #   combos: a vector of linear combinations, the subvectors have length
  #     J_i corresponding to regindex, and a first element which is Q
  #     BUT: if set to NULL, then all regressors in regindex are constrained
  #       to have the same value.
  #	Outputs:
  #		constraint: row vector of the form [Q , C], with C (constraint.mat)
  #     the vector of constraints and Q (constraint.vec) the  
  #     constraint constant, such that C psi = Q.
  #     When combos=NULL, constraint will be a matrix
  #	Requires: sigex.default, sigex.par2psi
  #
  ####################################################################
  
  par.mle <- sigex.default(mdl,data.ts,NULL)
  psi.mle <- sigex.par2psi(par.mle,mdl)
  psi.index <- NULL
  mark <- length(psi.mle) - length(par.mle[[4]]) 
  for(i in 1:N)
  {
    psi.index <- c(psi.index,mark + regindex[[i]])
    mark <- mark + dim(mdl[[4]][[i]])[2]
  }
  if(length(combos)>0) 
  {
    constraint <- rep(0,length(psi.mle))
    constraint[psi.index] <- combos[-1]
    constraint <- c(combos[1],constraint)
  } else
  {
    psi.one <- psi.index[1]
    constraint <- NULL
    for(k in 2:length(psi.index))
    {
      constraint.new <- rep(0,length(psi.mle))
      constraint.new[psi.one] <- 1
      constraint.new[psi.index[k]] <- -1 
      constraint.new <- c(0,constraint.new)
      constraint <- rbind(constraint,constraint.new)
    }
  }
  
  return(constraint)
}
