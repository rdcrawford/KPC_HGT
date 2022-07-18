#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
std::string MaskRecPositions( 
  std::string seq, std::vector< int > start, std::vector< int > end 
  ) 
{
  for ( int j = 0; j < start.size(); j++ )
  {
    for ( int i = start[j] - 1; i < end[j] - 1; i++ )
    {
      seq[i] = '-';
    }
  }
  return seq;
}

