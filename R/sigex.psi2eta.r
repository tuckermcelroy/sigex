sigex.psi2eta <- function(psi,constraint)
{
  
  ##########################################################################
  #
  #	sigex.psi2eta
  # 	    Copyright (C) 2019  Tucker McElroy
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
  #	Purpose: transform psi to eta
  #	Background:	
  #		psi refers to a vector of real numbers containing all
  #		hyper-parameters (i.e., reals mapped bijectively to the parameter
  #		manifold), whereas eta consists of free variables of real numbers
  #   that are mapped to a fixed set nu, and psi is a permuted
  #   vector of eta and nu.
  #	Notes: this is a functional inverse to sigex.eta2psi
  #	Inputs:
  #		psi: see background.  Must have length given by 
  #     dim(constraint.mat)[2] 
  #		constraint: matrix of the form [Q , C], with C (constraint.mat)
  #     the matrix of constraints and Q (constraint.vec) the vector
  #     of constraint constants, such that C psi = Q.
  #	Outputs:
  #		list of nu and eta: see background.
  #
  ####################################################################
  
  nu <- NULL
  eta <- psi
  if(length(constraint) > 0) 
  {
    constraint.mat <- constraint[,-1,drop=FALSE]
    constraint.vec <- constraint[,1,drop=FALSE]
    
    ## compute decomposition
    constraint.qr <- qr(constraint.mat)
    constraint.q <- qr.Q(constraint.qr)
    constraint.r <- qr.R(constraint.qr)
    constraint.pivot <- constraint.qr$pivot
    constraint.ipivot <- sort.list(constraint.pivot)
    
    ## compute mapping from free variables
    fixed.dim <- dim(constraint.mat)[1]
    free.dim <- dim(constraint.mat)[2] - fixed.dim 
    psi <- psi[constraint.pivot]
    nu <- psi[1:fixed.dim]
    eta <- psi[(fixed.dim+1):(fixed.dim+free.dim)]
   }
  
  return(list(nu,eta))
}


