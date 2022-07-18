#include <Rcpp.h>
#include <fstream>
using namespace std;

// [[Rcpp::plugins(cpp11)]]

// -----------------------------------------------------------------------------
// Get Reverse compliment
// Ryan D. Crawford
// 02/09/2019
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
std::string GetReverseCompliment( string &sequence )
{
  // Reverse the iput string
  reverse( sequence.begin(), sequence.end() );

  // Make sure the string is uppercase
  if (isupper(sequence[0]))
  {
    for (int i = 0; i < sequence.length(); i ++)
    {
      if (sequence[i] == 'A') sequence[i] = 'T';
      else if (sequence[i] == 'T') sequence[i] = 'A';
      else if (sequence[i] == 'C') sequence[i] = 'G';
      else if (sequence[i] == 'G') sequence[i] = 'C';
      else sequence[i] = 'N';
    }
  }
  else // transform the string in lower case
  {
    for (int i = 0; i < sequence.length(); i ++)
    {
      if (sequence[i] == 'a') sequence[i] = 't';
      else if (sequence[i] == 't') sequence[i] = 'a';
      else if (sequence[i] == 'c') sequence[i] = 'g';
      else if (sequence[i] == 'g') sequence[i] = 'c';
      else sequence[i] = 'n';
    }
  }
  return sequence;
}

// -----------------------------------------------------------------------------
                                                                                                     