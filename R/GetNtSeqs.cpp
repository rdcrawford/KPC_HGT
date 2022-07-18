#include <Rcpp.h>
#include <fstream>
using namespace std;

// [[Rcpp::plugins(cpp11)]]

// -----------------------------------------------------------------------------
// Parse Fasta Files By Position
// Ryan D. Crawford
// 02/09/2019
// -----------------------------------------------------------------------------
// This function takes as input a fasta file and vectors with start and end
// positions and writes the files to an output file stream.
// -----------------------------------------------------------------------------

// ---- Function Prototypes ----------------------------------------------------

// An nt sequence is input (passed by reference) and the sequence is updated
// at each position to the reverse complement
void GetReverseCompliment( string &sequence );

// ---- Define Function to Export ----------------------------------------------

// [[Rcpp::export]]
std::vector< std::string > GetNtSeqs(
  const Rcpp::DataFrame  &gffData,
  const vector< string > &genomeSeq
  )
{
  // ---- Constant Declarations ------------------------------------------------

  // const int ANNOT_ID  = 0; // Col index of annot IDs
  // const int GENE_NAME = 1; // Col index of the gene names
  const int CONTIG    = 2; // Col index of the contig index
  const int LEFT_POS  = 3; // Col index of Gene start positons
  const int RIGHT_POS = 4; // Col index of Gene end positions
  const int STRAND    = 5; // Col index indicating is on the antisense strand

  // ---- Variable Initializations ---------------------------------------------

  // Extract the vectors we need from the dataframe
  vector< string > ntSeqs( gffData.nrow() );         // Vector to store genes
  vector< int >    startPos  = gffData[ LEFT_POS ];  // Gene start positions
  vector< int >    endPos    = gffData[ RIGHT_POS ]; // Gene end positions
  vector< int >    contig    = gffData[ CONTIG ];    // Indexes of the contig
  vector< string > strand    = gffData[ STRAND ];    // Gene on +/- strand
  string           seq;                              // Gene coding sequence
  int              geneLen;                          // Length of the gene

  // ---- Parse the Fasta Files ------------------------------------------------

  // Open the output file stream
  assert( gffData.nrow( ) != 0 );

  // For each gene, write the sequence to the output file
  for ( int i = 0; i < gffData.nrow(); i++ )
  {
    geneLen = endPos[i] - startPos[i] + 1;

    // Subset the contig to extract the fragment
    seq = genomeSeq[ contig[i] ].substr( startPos[i], geneLen );

    // If this is the reverse strand, get the reverse compliment
    if ( strand[i] == "-" ) GetReverseCompliment( seq );

    ntSeqs[i] = seq;
  }

  return ntSeqs;
}

// -----------------------------------------------------------------------------

void GetReverseCompliment( string &sequence )
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
}

// -----------------------------------------------------------------------------
