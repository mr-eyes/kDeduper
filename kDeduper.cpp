#include "colored_kDataFrame.hpp"
#include <string>
#include <iostream>
#include <vector>
#include "kseqReader.hpp"
#include <sys/stat.h>
#include <map>
#include <fstream>
#include <cassert>
#include "Utils/kmer.h"

using namespace std;

inline bool file_exists(const std::string &name) {
    struct stat buffer
            {
            };
    return (stat(name.c_str(), &buffer) == 0);
}

uint64_t get_representativ_kmer(string &kmer_1, string &kmer_2) {
    uint64_t canonical_kmer_1 = kmer::str_to_canonical_int(kmer_1);
    uint64_t canonical_kmer_2 = kmer::str_to_canonical_int(kmer_2);
    if (canonical_kmer_1 < canonical_kmer_2) return canonical_kmer_1;
    else return canonical_kmer_2;
}

void debug_print_kmersCount(flat_hash_map<uint64_t, uint32_t> & kmersCount){
    cerr << "----- DEBUG START ------" << endl;
    for(const auto & item : kmersCount){
        cerr << item.first << ": " << item.second << endl;
    }
    cerr << "-----  DEBUG END  ------" << endl;
}

int main(int argc, char **argv) {

//    if (argc != 3) {
//        cerr << "run ./kDeduper <R1> <R2>" << endl;
//        exit(1);
//    }

    string R1_file = "/home/mabuelanin/Desktop/dev-plan/kDeduper/sample_data/R1.fastq"; // argv[1];
    string R2_file = "/home/mabuelanin/Desktop/dev-plan/kDeduper/sample_data/R2.fastq"; // argv[2];

    if (!file_exists(R1_file)) {
        throw std::runtime_error("Could not open R1 file");
    }

    if (!file_exists(R2_file)) {
        throw std::runtime_error("Could not open R2 file");
    }

    int kSize = 31;
    cerr << "counting number of reads ..." << endl;
    int count = 0;
    string line;
    ifstream file(R1_file);
    while (getline(file, line)) count++;
    int no_of_sequences = count / 2;

    // First iteration

    gzFile fp_1, fp_2;
    kseq_t *kseq_1, *kseq_2;
    fp_1 = gzopen(R1_file.c_str(), "r");
    fp_2 = gzopen(R2_file.c_str(), "r");

    kseq_1 = kseq_init(fp_1);
    kseq_2 = kseq_init(fp_2);

    flat_hash_map<uint64_t , uint32_t> kmers_to_count;

    int terminal_kmer_offset = 2;

    for (int seqCounter = 0; kseq_read(kseq_1) >= 0 && kseq_read(kseq_2) >= 0; seqCounter++) {

        uint32_t seq_1_length = string(kseq_1->seq.s).size();
        uint32_t seq_2_length = string(kseq_2->seq.s).size();

        if (seq_1_length < kSize || seq_2_length < kSize) continue;

        // Extract the start and end kmers
        string start_kmer = string(kseq_1->seq.s).substr(0 + terminal_kmer_offset, kSize);
        string end_kmer = string(kseq_2->seq.s).substr(seq_2_length - kSize - terminal_kmer_offset, kSize);

        // Get the representative kmer and increment its count
        uint64_t rep_kmer = get_representativ_kmer(start_kmer, end_kmer);
        kmers_to_count[rep_kmer]++;

    }

    debug_print_kmersCount(kmers_to_count);

}
