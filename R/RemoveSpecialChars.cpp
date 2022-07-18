#include <Rcpp.h>
using namespace Rcpp;


bool GetAlgnChar( std::vector< char > &algnChars, char inChar )
{
  auto it = find( algnChars.begin(), algnChars.end(), inChar );
  if ( it == algnChars.end() ) return true;
  return false;
}

// [[Rcpp::export]]
std::string RemoveSpecialChars( std::string seq )
{
  std::vector< char > algnChars = { 'A', 'C', 'G', 'T', 'N', '-' };
  for ( int i = 0; i < seq.size(); i++ )
  {
    if ( GetAlgnChar( algnChars, seq[i] ) ) seq[i] = 'N';
  }
  return seq;
}
