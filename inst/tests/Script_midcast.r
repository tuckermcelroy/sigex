# Scripts to do testing of mvar midcast rcpp code

library(Rcpp)

Rcpp::sourceCpp("/home/tucker-mcelroy/Documents/GitHub/sigex/inst/tests/mvar_midcast.cpp")

## Test values

z <- matrix(seq(1,15),nrow=3,ncol=5)
z[2,1] <- 1i
z[c(1,2),3] <- 1i

x.acf <- array(0,c(3,5,3))
x.acf[,1,] <- 2*diag(3) + c(1,2,-1) %x% t(c(1,2,-1))

delta <- c(1,-2,1)
delta <- c(1,-1)
#delta <- c(1,0,0,-1)

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

// DEPRECATE
// IntegerVector ind_data = cast_indices;
// ind_data.push_front(0);
// ind_data.push_back(T+1);
// IntegerVector gaps = diff(ind_data);
// gaps.push_back(0);
// if(max(gaps) > d) 
// { 
//  t_hash = min(subsetting(ind_data,gaps,d));
// } else
// {
//  Rcerr << "There are not " << d << " contiguous values\n";
//}  
