#ifndef BIO_H
#define BIO_H

#include <string>
#include <vector>


/*
This function should return true if and only if
every character in the input is one of ATCG.
*/
bool IsValidDNASequence(const std::string & input);

/*
This function should calculate the reverse complement DNA sequence.

The first argument is the sequence, the second argument is a pointer to
an empty string, which you should modify to store the result.

This is obtained by reversing the input sequence and swaping each
nucleotide/letter with it's complement:
A <-> T
C <-> G

Example:
input = AAATTCGGGG
reverse = GGGGCTTAAA
reverse complement = CCCCGAATTT
*/
void GetReverseComplementSequence(const std::string & input,  std::string * const output);


/*
This function should return the RNA transcript from a DNA sequence.

A RNA transcript is the reverse complement of the DNA sequence, but RNA
has U (uracil) instead of T (thiamine).

Make sure you don't have any redundant code.
*/
std::string GetRNATranscript(const std::string & input);


/*
This function should return a vector of vector of strings with each possible RNA
reading frame from the given DNA sequence.

There are three possible reading frames (because the genetic code has three
nucleotides per amino acid) in each direction (you can also transcribe DNA in
the reverse complement direction, called the antiparallel strand).

Order the sequences like so:
1: Original (0 offset)
2: Original (1 offset)
3: Original (2 offset)
4: Antiparallel (0 offset)
5: Antiparallel (1 offset)
6: Antiparallel (2 offset)

With in the input sequence of: AATTCCCGAAA
Original RNA transcript = AAUUCCCGAAA
Antiparallel RNA transcript = UUUGCCCAAUU

The offsets (starting at pos 0, 1, and 2) of the two RNA transcripts
UUUGCCCAAUU
UUGCCCAAUU
UGCCCAAUU
AAUUCCCGAAA
AUUCCCGAAA
UUCCCGAAA

Instead of returning a vector of 6 strings, break each string into a vector
of length 3 strings (called codons) These codons will be useful
for the next translation step.

UUUGCCCAAUU -> {"UUU", "GCC", "CAA"}
// drop any remaining letters that don't fill a codon
UUGCCCAAUU -> {"UUG", "CCC", "AAU"}
UGCCCAAUU -> {"UGC", "CCA", "AUU"}
AAUUCCCGAAA -> {"AAU", "UCC", "CGA"}
AUUCCCGAAA -> {"AUU", "CCC", "GAA"}
UUCCCGAAA -> {"UUC", "CCG", "AAA"}
*/
std::vector< std::vector<std::string> > GetReadingFramesAsCodons(const std::string & input);

/*
This function translates/converts a vector<string> (vector of codons) into a
string of amino acids using the genetic code
(see https://en.wikipedia.org/wiki/Genetic_code).

For example, the codons:
{"UUU", "GCC", "CAA"}
translates to:
F (Phenylalanine), A (Alanine), Q (Glutamine)
abreviated:
FAQ

http://www.soc-bdr.org/rds/authors/unit_tables_conversions_and_genetic_dictionaries/genetic_code_tables/


To make your lives easier, here's a list of the possible codons:
"GCU", "GCC", "GCA", "GCG", "CGU", "CGC", "CGA", "CGG", "AGA", "AGG",
"AAU", "AAC", "GAU", "GAC", "UGU", "UGC", "CAA", "CAG", "GAA", "GAG",
"GGU", "GGC", "GGA", "GGG", "CAU", "CAC", "AUU", "AUC", "AUA", "UUA",
"UUG", "CUU", "CUC", "CUA", "CUG", "AAA", "AAG", "AUG", "UUU", "UUC",
"CCU", "CCC", "CCA", "CCG", "UCU", "UCC", "UCA", "UCG", "AGU", "AGC",
"ACU", "ACC", "ACA", "ACG", "UGG", "UAU", "UAC", "GUU", "GUC", "GUA",
"GUG", "UAG", "UGA", "UAA"

And there corresponding amino acids ("*" represents STOP codons,
more on them later):

"A", "A", "A", "A", "R", "R", "R", "R", "R", "R", "N", "N", "D", "D",
"C", "C", "Q", "Q", "E", "E", "G", "G", "G", "G", "H", "H", "I", "I",
"I", "L", "L", "L", "L", "L", "L", "K", "K", "M", "F", "F", "P", "P",
"P", "P", "S", "S", "S", "S", "S", "S", "T", "T", "T", "T", "W", "Y",
"Y", "V", "V", "V", "V", "*", "*", "*"
*/
std::string Translate(const std::vector<std::string> & codon_sequence);

/*
This function takes a DNA sequence and returns the longest possible
amino acid sequence / protein that is encoded by that sequence
(open reading frame). A valid open reading frame begins with the
codon AUG (the amino acid, Methionine (M)) and runs until a stop codon (*)
is encountered. There may be multiple open reading frames in a sequence, and
you need to check all six reading frames in order given by
get_reading_frames_as_codons. If there are ties for longest, favor the first
one found.

Return the longest open reading frame as an amino acid sequence. It must start
with an 'M' and end with a '*' with no other '*''s within.
*/
std::string GetLongestOpenReadingFrame(const std::string & DNA_sequence);





#endif
