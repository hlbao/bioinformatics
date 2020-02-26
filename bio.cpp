#include "bio.h"
#include <iostream>
#include <string>
#include <vector>
using std::string;

//This function aims to determine whether input is ATGC, if yes return true, if no return false;
bool IsValidDNASequence(const std::string & input){
    for(int i = 0; i < static_cast<int>(input.size()); i++){
        if(input[i] != 'A' && input[i] != 'T' && input[i] != 'G' && input[i] != 'C'){
            return false;
        }
    }
    return true;
}

void GetReverseComplementSequence(const std::string & input,  std::string & const output){
    // iterates through string from end to beginging and adds
    // the complement of the current element to a new string
    for(int i = static_cast<int>(input.size()) -1; i>=0; i--){
        if(input[i] == 'A'){
            *output+='T';
        }
        else if(input[i] == 'T'){
            *output += 'A';
        }
        else if(input[i] == 'G'){
            *output+= 'C';
        }
        else{
            *output += 'G';
        }
    }
}

std::string GetRNATranscript(const std::string & input){
    // iterates through the input string from end to start
    // adds complement of element to a new string
    std::string RNA;
     for(int i = static_cast<int>(input.size()) -1; i>=0; i--){
        if(input[i] == 'A'){
            RNA +='U';
        }
        else if(input[i] == 'T'){
            RNA += 'A';
        }
        else if(input[i] == 'G'){
            RNA += 'C';
        }
        else{
            RNA += 'G';
        }
    }
    return RNA;
}

std::vector<std::vector<std::string>> GetReadingFramesAsCodons(const std::string & input){
    std::string compDNA;
    GetReverseComplementSequence(input, &compDNA);

    std::string RNA = GetRNATranscript(input);
    std::string compRNA = GetRNATranscript(compDNA);

    std::vector<std::vector<std::string>> codons(6);
    // iterates through the (normal) RNA and creates sections of 3 and adds it to the vector of vectors
    for(int offset = 0; offset < 3; offset++){
        for(int i = offset; i + 2 < static_cast<int>(RNA.size()); i+=3){
            codons[offset].push_back(RNA.substr(i, 3));
        }
    }
    // iterates through the (complementary) RNA and creates sections of 3 and adds it to the vector of vectors
    for(int offset = 0; offset < 3; offset++){
        for(int i = offset; i + 2 < static_cast<int>(compRNA.size()); i+=3){
            codons[offset + 3].push_back(compRNA.substr(i, 3));
        }
    }
    return codons;
}

std::string Translate(const std::vector<std::string> & codon_sequence){
    //hard-coded codons and amino acids
    std::vector<std::string> codons = {"GCU", "GCC", "GCA", "GCG", "CGU", "CGC",
    "CGA", "CGG", "AGA", "AGG","AAU", "AAC", "GAU", "GAC", "UGU", "UGC",
    "CAA", "CAG", "GAA", "GAG", "GGU", "GGC", "GGA", "GGG", "CAU", "CAC",
    "AUU", "AUC", "AUA", "UUA","UUG", "CUU", "CUC", "CUA", "CUG", "AAA",
    "AAG", "AUG", "UUU", "UUC", "CCU", "CCC", "CCA", "CCG", "UCU", "UCC",
    "UCA", "UCG", "AGU", "AGC", "ACU", "ACC", "ACA", "ACG", "UGG", "UAU",
    "UAC", "GUU", "GUC", "GUA", "GUG", "UAG", "UGA", "UAA"};

    std::vector<std::string> amino_acids = {"A", "A", "A", "A", "R", "R",
    "R", "R", "R", "R", "N", "N", "D", "D", "C", "C", "Q", "Q", "E", "E",
    "G", "G", "G", "G", "H", "H", "I", "I", "I", "L", "L", "L", "L", "L",
    "L", "K", "K", "M", "F", "F", "P", "P", "P", "P", "S", "S", "S", "S",
    "S", "S", "T", "T", "T", "T", "W", "Y", "Y", "V", "V", "V", "V", "*",
    "*", "*"};

    std::string aa_sequence;
    // goes through the passed codon sequence and finds the index of it in the
    // hard-coded vector of codons. Uses the found index to add the respective
    // amino acid to the chain that will eventually be returned.
    for(std::string codon : codon_sequence){
        for(int i = 0; i < static_cast<int>(codons.size()); i++){
            if (codon == codons[i]){
                aa_sequence += amino_acids[i];
            }
        }
    }
    return aa_sequence;
}

std::string GetLongestOpenReadingFrame(const std::string & DNA_sequence){
    // first, create the codon sequences
    std::vector<std::vector<std::string>> codon_sequences = GetReadingFramesAsCodons(DNA_sequence);
    //initialize the longest sequence string to be used later
    std::string longest_sequence;
    for(int i = 0; i < static_cast<int>(codon_sequences.size()); i++){
        // iterates through the 6 codon sequences in the creaetd vector
        std::string aa_sequence = Translate(codon_sequences[i]); // trasnaltes current codon sequence into aimino acid sequence
        bool create = false;
        std::string current_sequence;
        for(int j = 0; j < static_cast<int>(aa_sequence.size()); j++){ // iterates through amino acid sequence
            if(aa_sequence[j] == 'M'){ // if an M is identified, start constructing the current chain
                current_sequence += 'M';
                create = true;
            }
            if(create && aa_sequence[j] != 'M'){ // constructs current chain if M was found
                current_sequence += aa_sequence[j];
                if(aa_sequence[j] == '*'){ // ends construction and checks if the built chain is longer than the stored longest chain
                    create = false;
                    if(static_cast<int>(current_sequence.size()) > static_cast<int>(longest_sequence.size())){
                        longest_sequence = current_sequence; // overrieds longest chain if new chain has a longer length
                    }
                    current_sequence = ""; // resets the current chain
                }
            }

        }
    }
    return longest_sequence;
}
