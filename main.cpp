#include <iostream>
using std::cout; using std::cin; using std::endl;
#include <string>
using std::string;
#include <sstream>
using std::istringstream;
#include <vector>
using std::vector;

#include "bio.h"

int main() {
    cout << "Project 2 (Genetics)" << endl;

    cout << "Please provide a DNA sequence:" << endl;
    string sequence;
    cin >> sequence;

    cout << "The DNA sequence you provided is:" << endl
         << sequence << endl;

    if (! IsValidDNASequence(sequence)) {
        cout << "The above sequence isn't a valid DNA sequence." << endl;
        return 1;
    }

    cout << "The reverse complement DNA sequence is:" << endl;
    string reverse_complement_sequence;
    GetReverseComplementSequence(sequence, &reverse_complement_sequence);
    cout << reverse_complement_sequence << endl;

    cout << "The RNA transcript for the original sequence is:" << endl;
    cout << GetRNATranscript(sequence) << endl;

    cout << "The codon reading frames are:" << endl;
    for (vector<string> frame : GetReadingFramesAsCodons(sequence)) {
        for (string codon : frame) {
            cout << codon << ", ";
        }
        cout << endl;
    }
    cout << "The amino acid sequence for each of the above are:" << endl;
    for (vector<string> frame : GetReadingFramesAsCodons(sequence)) {
        cout << Translate(frame) << endl;
    }

    string longest_open_reading_frame = GetLongestOpenReadingFrame(
        sequence);
    if (longest_open_reading_frame.empty()) {
        cout << "There is no open reading frame." << endl;
    } else {
        cout << "The longest open reading frame is:" << endl
             << longest_open_reading_frame << endl;
    }


}
