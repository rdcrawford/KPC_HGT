#include <Rcpp.h>
#include <fstream>
using namespace std;

// [[Rcpp::plugins(cpp11)]]

// -----------------------------------------------------------------------------
// Count gap positions
// Ryan D. Crawford
// 02/09/2019
// -----------------------------------------------------------------------------
// This function takes as input a fasta file and vectors with start and end
// positions and writes the files to an output file stream.
// -----------------------------------------------------------------------------

// ---- Define Function to Export ----------------------------------------------

// [[Rcpp::export]]
std::vector< int > CountGapPositions( const vector< string > &algn )
{
  const int          nSeqs = algn.size();
  std::vector< int > gapCounts( nSeqs );
  int                seqCount = 0;
  
  for ( int i = 0; i < nSeqs; i++ )
  {
    for ( auto it = algn[i].begin(); it != algn[i].end(); it++)
      if ( *it == '-' ) gapCounts[i] ++;
  }
  return gapCounts;
}

// -----------------------------------------------------------------------------