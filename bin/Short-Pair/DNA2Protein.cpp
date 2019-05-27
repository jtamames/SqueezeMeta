// DNA2Protein.cpp : Defines the entry point for the console application.
//
#include <cmath>
#include <cstdio>
#include <cassert>
#include <cstring>
#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <fstream>
#include <iostream>
using namespace std;
class DNASeq {
 public:
	DNASeq(const string& seq_name, const string& seq);
	string seq_name() const {return seq_name_;}
	string seq() const {return seq_;}
	int length() const {return length_;}
	string reverse_complement() const {return reverse_complement_;}
	vector<string> protein_seqs(const int& first_frame, const int& last_frame) const; // Return 6-frame translations of seq_
	friend int Codon2AminoAcidIndex(const char* codon);
	friend char Codon2AminoAcid(const char* codon);

 private:
  string seq_name_; 
	string seq_; 
	int length_;
	string reverse_complement_; //  The reverse complement of the input sequence.
	static const int kCodonNum_;
	static const char* kCodonTable_; // Amino acids for 64 possible codons plus 'X' for codons wtih unexpected characters other than A, C, G and T
	static const int kCodonIndexTable_[65]; // Indices for each amino acid according to their position in hmm files
	static const int kMaxLineLength_;
};
int Codon2AminoAcidIndex(const char* codon);
char Codon2AminoAcid(const char* codon);
void OpenOutFiles(const char* out_file_base_name, const int& first_frame, const int& last_frame, ofstream* out_files);
void OutputOneSeqTranslation(const DNASeq& dna, const int& first_frame, const int& last_frame, ofstream* out_files);
void CloseOutFiles(ofstream* out_files, int size);
int main(int argc, char* argv[]) {
	if (argc != 4) {
		cerr << "DNA2Protein [frame optionsuch as 1-3 or 4-6] [fasta file] [output fasta file base name]" << endl;
		exit(1);
	}
  if (strlen(argv[1]) != 3 || argv[1][1] != '-') {
    cerr << "The format of the frame option should be like 1-3 or 4-6 or 1-6" << endl;
	  exit(1);
  }
  int first_frame, last_frame; 
  first_frame = argv[1][0] - 48;
  last_frame = argv[1][2] - 48;
  int frame_size = last_frame - first_frame + 1;
	if (first_frame < 1 || last_frame >6 || first_frame > last_frame) {
    cerr << "Frames can only be from 1-6 and the smaller frame should come first" << endl;
		exit(1);         
	}
	string line;
	const char* fasta_file_name = argv[2];
	ifstream fasta_file(fasta_file_name);
	assert(fasta_file.is_open());
	string seq, seq_name;
	ofstream out_files[6];
  OpenOutFiles(argv[3], first_frame, last_frame, out_files);
	while (getline(fasta_file, line)) {
		if (line.empty()) {
			continue;
		}
		if (line[0] == '>') {
			if (!seq_name.empty()) {
				DNASeq dna(seq_name, seq);
				OutputOneSeqTranslation(dna, first_frame, last_frame, out_files);
			}
			seq_name = line.substr(1);
			seq.clear();
		} else {
			seq.append(line);
		}
	}
	DNASeq dna(seq_name, seq);
	OutputOneSeqTranslation(dna, first_frame, last_frame, out_files);
	CloseOutFiles(out_files, frame_size);
	fasta_file.close();
	return 0;
}
const int DNASeq::kMaxLineLength_ = 2048;
const int DNASeq::kCodonNum_ = 64;
const char* DNASeq::kCodonTable_ = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVXYXYSSSSXCWCLFLF";
const int DNASeq::kCodonIndexTable_[65] = {8,11,8,11,16,16,16,16,14,15,14,15,7,7,10,7,13,6,13,6,12,12,12,12,14,14,14,14,9,9,9,9,3,2,3,2,0,0,0,0,5,5,5,5,17,17,17,17,20,19,20,19,15,15,15,15,20,1,18,1,9,4,9,4,20};
DNASeq::DNASeq(const string& seq_name, const string& seq)
    : seq_(seq),
      seq_name_(seq_name)
{
	seq_name_ = seq_name;
	length_ = seq_.length();
	reverse_complement_ = seq_;
	for (int i = 0; i < length_; i++) {
		if (seq_[i] > 'Z') {		
			seq_[i] = seq_[i] - 32;
		}
		switch (seq_[i]) {
			case 'A':
				reverse_complement_[length_ - 1 - i]= 'T';
				break;
			case 'C':
				reverse_complement_[length_ - 1 - i] = 'G';
				break;
			case 'G':
				reverse_complement_[length_ - 1 - i] = 'C';
				break;
			case 'T':
				reverse_complement_[length_ - 1 - i] = 'A';
				break;
			default:
				reverse_complement_[length_ - 1 - i] = seq_[i];
				break;
		}
	}
}
//frame translation
vector<string> DNASeq::protein_seqs(const int& first_frame, const int& last_frame) const {
	vector<string> protein_seqs(last_frame - first_frame + 1);
	const char* seq = seq_.c_str();
  const char* reverse_complement = reverse_complement_.c_str();
  const char* translate_seq;
	for (int i = first_frame - 1; i < last_frame; ++i) {
    int frame_index = i + 1 - first_frame; //  index of the current frame in the vector
		int j = i;
    if (i < 3) {
      translate_seq = seq;
    } else {
			translate_seq = reverse_complement;
			j -= 3;
		}	
		for (; j + 3 <= length_; j += 3) {
			protein_seqs[frame_index].push_back(Codon2AminoAcid(translate_seq + j));
		}   
	}
	return protein_seqs;
}
int Codon2AminoAcidIndex(const char* codon) { // Convert a codon into its corresponding amino acid character
		int index = 0;
		for (int i =0; i < 3; i++) {
			switch (codon[2 - i]) {
				case 'A':
					//index += pow(4.,i) * 0;
					break;
				case 'C':
					index += (1<<(i<<1));
					break;
				case 'G':
					//index += pow(4.,i) * 2;
					index |= (1<<(i<<1))<<1;
					break;
				case 'T': 
					//index += pow(4.,i) * 3;
					index |= (1<<((i<<1)|1))+(1<<(i<<1));
					break;
				default:
					return 20; // Invalid characters lead to X amino acid
			}
		}
		return DNASeq::kCodonIndexTable_[index];
}
char Codon2AminoAcid(const char* codon) { // Convert a codon into its corresponding amino acid character
		double index = 0;
		for (int i =0; i < 3; i++) {
			switch (codon[2 - i]) {
				case 'A':
					index += pow(4.,i) * 0;
					break;
				case 'C':
					index += pow(4.,i) * 1;
					break;
				case 'G':
					index += pow(4.,i) * 2;
					break;
				case 'T': 
					index += pow(4.,i) * 3;
					break;
				default:
					return 'X'; // Invalid characters lead to X amino acid
			}
		}
		return DNASeq::kCodonTable_[static_cast<int>(index)];
	}
void OpenOutFiles(const char* out_file_base_name, const int& first_frame, const int& last_frame, ofstream* out_files) {
  char out_file_name[100];
 	char suffix_name[8] = ".framex";
  int first_index = 0; // index of the first frame in the vector, always begins with 0
  int last_index = last_frame - first_frame; // index of the last frame in the vector
  for (int i = first_index; i <= last_index; ++i) {
    strcpy(out_file_name, out_file_base_name);
		suffix_name[6] = i + first_frame + 48; // the file name should be the same with the frame instead of index
    strcat(out_file_name, suffix_name);
    out_files[i].open(out_file_name);
    assert(out_files[i].is_open());
	}
}
void OutputOneSeqTranslation(const DNASeq& dna, const int& first_frame, const int& last_frame, ofstream* out_files) {
	vector<string> proteins = dna.protein_seqs(first_frame, last_frame);
  int size = last_frame - first_frame + 1;
	for (int i = 0; i < size; ++i) {
 		out_files[i] << ">" << dna.seq_name() << "\n" << proteins[i] << endl;
  }
}
void CloseOutFiles(ofstream* out_files, int size) {
	for (int i = 0; i < size; ++i) {
	  	out_files[i].close();
	}
}

