#include <RcppArmadillo.h>
#include <progress.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif


using namespace std;
using namespace Rcpp;


// [[Rcpp::export]]
arma::sp_mat cellDist(const arma::sp_mat& m, int ncores=1,bool verbose=true, bool full=false, bool diag=true,bool absolute=false)
{
  typedef arma::sp_mat::const_col_iterator iter;
  arma::sp_mat d(m.n_cols,m.n_cols);
  int tot = (full) ? m.n_cols*m.n_cols - ((diag) ? m.n_cols : 0) : (m.n_cols*m.n_cols + ((diag) ? m.n_cols : 0))/2;
  Progress p(tot, verbose);
#pragma omp parallel for num_threads(ncores) shared(d)
  for(unsigned int i=0;i<m.n_cols;i++) {
    for(unsigned int j = diag ? i : i+1;j<m.n_cols;j++) {
      if ( !Progress::check_abort())
      {
        p.increment(); // update progress
        iter i_iter = m.begin_col(i);
        iter j_iter = m.begin_col(j);
        double dem=0;

        while( (i_iter != m.end_col(i)) && (j_iter != m.end_col(j)) )
        {
          if(i_iter.row() == j_iter.row())
          {
            dem+= (absolute) ? std::abs((*i_iter)+(*j_iter)) : std::abs(arma::sign(*i_iter)+arma::sign(*j_iter));
            ++i_iter;
            ++j_iter;
          } else {
            if(i_iter.row() < j_iter.row())
            {
              dem+=std::abs((*i_iter));
              ++i_iter;
            } else {
              dem+=std::abs((*j_iter));
              ++j_iter;
            }
          }
        }
        for(; i_iter != m.end_col(i); ++i_iter) { dem+=std::abs(*i_iter); }
        for(; j_iter != m.end_col(j); ++j_iter){ dem+=std::abs(*j_iter); }
        d(j,i) = (2*m.n_cols - dem)/(2*m.n_cols);
        if (full) {d(i,j) = d(j,i);}
      }
    }
  }

  return(d);
}

