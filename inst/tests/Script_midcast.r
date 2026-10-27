# Scripts to do testing of mvar midcast rcpp code

library(Rcpp)

Rcpp::sourceCpp("/home/tucker-mcelroy/Documents/GitHub/sigex/inst/tests/mvar_midcast.cpp")

## Test values

z <- matrix(seq(1,15),nrow=3,ncol=5)
z[2,1] <- 1i
z[c(1,2),4] <- 1i

x.acf <- array(0,c(3,5,3))
x.acf[,1,] <- diag(3)

delta <- c(1,-1)


## Tests of helper functions

print(getRe(z[,1]))
print(getIm(z[,1]))
print(get_rag_index(z[,1]))

print(matmult(matrix(seq(1,4),nrow=2),matrix(seq(1,4),nrow=2)))
print(kronprod(matrix(seq(1,4),nrow=2),matrix(seq(1,4),nrow=2)))
print(vec2mat(seq(1,4)))

print(mvar_midcast(x.acf,z,delta))


#####################################################

## Pieces of main midcast rcpp function, stored here until compilation works

## NOTE: mvar_midcast should be type List eventually i.e.

//List mvar_midcast(arma::cube x_acf, ComplexMatrix z, NumericVector delta, bool debug = false) 

## Other code chunks:


if(max(gaps) > d) 
{ 
  t_hash = min(subsetting(ind_data,gaps,d));
} else
{
  Rcerr << "There are not " << d << " contiguous values\n";
}  

//  Rcout << "t.hash  "  << t_hash;

//  Main code block for nonstationary case (d > 0)      
if(d > 0) 
{
  
  // t = t_hash + d case as initialization
  //  get predictors based on observations t_hash+1 to t
  arma::mat l_pred (N,N*d);
  arma::vec del_aft = -1*rev(tail(delta,d))/delta(0);
  arma::mat del_aft_mat = vec2mat(del_aft);
  l_pred = kronprod(del_aft_mat,idmat);
  l_pred = l_pred.t();
  arma::mat l_derp (N,N*d);
  arma::vec del_fore = -1*rev(head(delta,d))/delta(d);
  arma::mat del_fore_mat = vec2mat(del_fore);
  l_derp = kronprod(del_fore_mat,idmat);
  l_derp = l_derp.t();
  arma::mat v_pred = x_acf.slice(0)/(delta(0)*delta(0));
  arma::mat v_derp = x_acf.slice(0)/(delta(d)*delta(d));
  
  //  get casts and covars of observations t_hash+1 to t based
  //    on sigma-field_{t_hash+1:t}
  //   Note: store more than needed in preds_x, makes it easier 
  //    for indexing later
  arma::mat preds_x (N,t_hash+d);
  for(int t = 0; t < t_hash+d; t++)
  {
    preds_x.col(t) = getRe(z(_,t));
  }
  arma::mat new_covar;
  arma::mat casts_x;
  arma::mat casts_var;
  arma::mat eps;
  double Qseq = 0;
  double logdet = 0;
  
  // track indices of casted variables up to present time
  IntegerVector init_index = seq(t_hash+1,t_hash+d);
  IntegerVector cast_index_t = intersect(cast_indices,init_index);
  int t_star = T;
  int t_len = d;
  
  // Forward Pass:	
    if(t_hash < T-d) 
    {
      for(int t = t_hash+d+1; t <= T; t++)
      { 
        int tm1 = t-1;
        // determine whether full info, or partial/completely missing
        //  base case: full info
        arma::mat select_mat;
        arma::mat omit_mat;
        LogicalVector rag_vec = ragged(_,tm1);          
        IntegerVector raggeds;
        IntegerVector non_raggeds;
        for(int k = 0; k < N; k++)
        {
          if(rag_vec(k)) { non_raggeds.push_back(k+1); } else { raggeds.push_back(k+1); }
        }  
        //        if(raggeds.length() > 0) { omit_mat = idmat.submat(raggeds,seq(1,N)); }
        //        if(non_raggeds.length() > 0) { select_mat = idmat.submat(non_raggeds,seq(1,N)); }
        //if(raggeds.length() > 0) { omit_mat = idmat.submat(1,2,1,1); }
        
        
        Rcout << "select matrix "  << select_mat << "\n";
        Rcout << "omit matrix "  << omit_mat << "\n";
        Rcout << "raggeds " << raggeds << "\n";       
        Rcout << "non raggeds " << non_raggeds << "\n";       
        
        
      } 
    }  // end forward pass clause
  
}   // end d > 0  clause

List L = List::create(ragged);

return L;

}





######################################################

## Deprecated Rcpp functions...just storing them here!

// DEPRECATE - unneeded?
  // // [[Rccp::export]]
// IntegerVector get_full_indices(ComplexMatrix z, IntegerVector all_indices)
// {
  //   int num_rows = z.nrow(); 
  //   int num_cols = z.ncol(); 
  //   
    //   // Vector to collect indices that meet the criteria
  //   std::vector<int> kept_indices;
  //   kept_indices.reserve(num_rows); // Pre-allocate maximum possible size for speed
  //   
    //   // Loop over each row (matches colSums behavior)
  //   for (int i = 0; i < num_rows; i++) {
    //     bool has_im_equal_to_one = false;
    //     ComplexVector zvec = z.row(i);
    //     arma::vec imPart = getIm(zval);
    //     
      //     // Loop over columns of that row
    //     for (int j = 0; j < num_cols; j++) {
      //       if(imPart(j) == 1)  // This means zval[j] is a missing value
      //       {
        //         has_im_equal_to_one = true;
        //         break; // Performance boost: Stop checking this column immediately!
          //       }
      //     }
    //     
      //     // If NO elements had an imaginary part of 1, keep this column's index
//     // Note: C++ uses 0-based indexing; ensure all_indices matches your intent
//     if (!has_im_equal_to_one) {
//       kept_indices.push_back(all_indices[j]);
//     }
//   }
//   
//   // Wrap back into an R-friendly IntegerVector to return to R
//   return wrap(kept_indices);
// }

// DEPRECATE - unneeded?
// // [[Rcpp::export]]
// IntegerVector get_cast_indices(IntegerVector all_indices, IntegerVector full_indices) {
//   // C++ set_difference requires inputs to be sorted
//   IntegerVector all_sorted = clone(all_indices);
//   IntegerVector full_sorted = clone(full_indices);
//   std::sort(all_sorted.begin(), all_sorted.end());
//   std::sort(full_sorted.begin(), full_sorted.end());
//   
//   // Create a destination vector container
//   std::vector<int> result;
//   
//   // Perform the set difference
//   std::set_difference(
//     all_sorted.begin(), all_sorted.end(),
//     full_sorted.begin(), full_sorted.end(),
//     std::back_inserter(result)
//   );
//   
//   // Convert back to Rcpp IntegerVector and return
//   return wrap(result);
// }
