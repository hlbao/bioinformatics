#include <string>
#include <map>
#include <vector>
#include <iostream>
using namespace std;

bool IsValidDNASequence(const std::string & input){
	for (int i =0 ; i < (int)input.size(); i++)
		if(input[i] != 'A' && input[i] != 'T' && input[i] != 'G' && input[i] != 'C')
			return 0;
	return 1;
}


void GetReverseComplementSequence(const string & input, string * const output){
	map<char, char> mp;
	mp['A'] = 'T';
	mp['T'] = 'A';
	mp['G'] = 'C';
	mp['C'] = 'G';
	for (int i = (int)input.size()-1; i >= 0; i--)
		*output += mp[input[i]];
}

string GetRNATranscript(const  string & input){
	map<char, char> mp;
	mp['A'] = 'U';
	mp['T'] = 'A';
	mp['G'] = 'C';
	mp['C'] = 'G';
	string output;
	for (int i = (int)input.size()-1; i >= 0; i--)
        output += mp[input[i]];
	return output;
}

vector<vector<string> > GetReadingFramesAsCodons(const  string & input){
    string cDNA;
    GetReverseComplementSequence(input, &cDNA);

    string RNA = GetRNATranscript(input);
    string cRNA = GetRNATranscript(cDNA);

    vector<vector<string> > codons(6);
    for(int offset = 0; offset < 3; offset++){
        for(int i = offset; i + 2 < (int)RNA.size(); i+=3){
            codons[offset].push_back(RNA.substr(i, 3));
        }
        for(int i = offset; i + 2 < (int)cRNA.size(); i+=3){
            codons[offset + 3].push_back(cRNA.substr(i, 3));
        }
    }
    // cout <<"!!!!!!!!!!\n";
    // cout << cRNA << endl;
    // string t = cRNA.substr(0, 3);
    // cout << int(t[0]) <<" " << int(t[1]) << " "<<int(t[2]) << endl;
    // cout << cRNA.substr(0, 3).size() << endl;
    // cout <<"!!!!!!!!!!\n";
    // for (int i =0 ; i < codons[3].size(); i++) cout << codons[3][i] <<" "; cout <<endl;
    return codons;
}

string Translate(const vector<string> & codon_sequence){
 	vector<string> codons = {"GCU", "GCC", "GCA", "GCG", "CGU", "CGC", 
    "CGA", "CGG", "AGA", "AGG","AAU", "AAC", "GAU", "GAC", "UGU", "UGC", 
    "CAA", "CAG", "GAA", "GAG", "GGU", "GGC", "GGA", "GGG", "CAU", "CAC", 
    "AUU", "AUC", "AUA", "UUA","UUG", "CUU", "CUC", "CUA", "CUG", "AAA", 
    "AAG", "AUG", "UUU", "UUC", "CCU", "CCC", "CCA", "CCG", "UCU", "UCC",
    "UCA", "UCG", "AGU", "AGC", "ACU", "ACC", "ACA", "ACG", "UGG", "UAU",
    "UAC", "GUU", "GUC", "GUA", "GUG", "UAG", "UGA", "UAA"};

    vector<string> amino_acids = {"A", "A", "A", "A", "R", "R",
    "R", "R", "R", "R", "N", "N", "D", "D", "C", "C", "Q", "Q", "E", "E",
    "G", "G", "G", "G", "H", "H", "I", "I", "I", "L", "L", "L", "L", "L",
    "L", "K", "K", "M", "F", "F", "P", "P", "P", "P", "S", "S", "S", "S",
    "S", "S", "T", "T", "T", "T", "W", "Y", "Y", "V", "V", "V", "V", "*",
    "*", "*"};

    map<string, string> mp;
    for (int i=0;i < 64; i++) mp[codons[i]] = amino_acids[i];

    string ans;
    for(string codon : codon_sequence)
    	ans += mp[codon];
    return ans;
 }

string GetLongestOpenReadingFrame(const  string & DNA_sequence){
 	vector<vector<string> > seq = GetReadingFramesAsCodons(DNA_sequence);
 	string ans;
 	for (vector<string> condons : seq){
 		string w = Translate(condons);        
        for (int i = 0, j;i < (int)w.size(); i ++)
            if (w[i] == 'M'){
                for (j = i + 1; j < (int)w.size(); j ++)
                    if(w[j] == '*'){
                        string tmp = w.substr(i, j - i + 1);
                        // cout << "!!!!" << tmp << endl;
                        if(tmp.size() > ans.size()) ans = tmp;
                        break;
                    }
                i = j;
            }
 		
 	}
 	return ans;
}

