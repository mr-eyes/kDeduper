#include <string>
#include <iostream>
#include <vector>
#include "kseqReader.hpp"
#include <sys/stat.h>
#include <map>
#include "kmer.h"
#include "phmap.h"
#include "fstream"
using namespace std;
using namespace phmap;

inline bool file_exists(const std::string &name) {
    struct stat buffer{};
    return (stat(name.c_str(), &buffer) == 0);
}

uint64_t get_representative_kmer(string &kmer_1, string &kmer_2) {
    uint64_t canonical_kmer_1 = kmer::str_to_canonical_int(kmer_1);
    uint64_t canonical_kmer_2 = kmer::str_to_canonical_int(kmer_2);
    if (canonical_kmer_1 < canonical_kmer_2) return canonical_kmer_1;
    else return canonical_kmer_2;
}

void debug_print_kmersCount(flat_hash_map<uint64_t, uint32_t> &kmersCount) {
    cerr << "----- DEBUG START (kmer count) ------" << endl;
    for (const auto &item : kmersCount) {
        cerr << item.first << ": " << item.second << endl;
    }
    cerr << "-----  DEBUG END  ------" << endl;
}

std::string reverse_complement(std::string seq) {
    auto lambda = [](const char c) {
        switch (c) {
            case 'A':
                return 'T';
            case 'G':
                return 'C';
            case 'C':
                return 'G';
            case 'T':
                return 'A';
            case 'N':
                return 'N';
            default:
                throw std::domain_error("Invalid nucleotide.");
        }
    };
    std::transform(seq.cbegin(), seq.cend(), seq.begin(), lambda);
    reverse(seq.begin(), seq.end());
    return seq;
}

class fastqWriter {

public:
    ofstream fileStream_r1, fileStream_r2;

    explicit fastqWriter(const string &filename_prefix) {
        this->fileStream_r1.open(filename_prefix + "_R1.fastq");
        this->fileStream_r2.open(filename_prefix + "_R2.fastq");
    }

    void write(kseq_t *seq1, kseq_t *seq2) {

        this->fileStream_r1 << '@';
        this->fileStream_r1 << seq1->name.s;
        this->fileStream_r1 << ' ';
        this->fileStream_r1 << seq1->comment.s;
        this->fileStream_r1 << endl;
        this->fileStream_r1 << seq1->seq.s;
        this->fileStream_r1 << "\n+\n";
        this->fileStream_r1 << seq1->qual.s;
        this->fileStream_r1 << endl;

        this->fileStream_r2 << '@';
        this->fileStream_r2 << seq2->name.s;
        this->fileStream_r2 << ' ';
        this->fileStream_r2 << seq2->comment.s;
        this->fileStream_r2 << endl;
        this->fileStream_r2 << seq2->seq.s;
        this->fileStream_r2 << "\n+\n";
        this->fileStream_r2 << seq2->qual.s;
        this->fileStream_r2 << endl;
    }

    void close() {
        fileStream_r1.close();
        fileStream_r2.close();
    }

};


uint64_t canonicalFragmentHash(kseq_t *kseq_1, kseq_t *kseq_2, hash<string> &fragmentHasher) {

    string fragment_seq = kseq_1->seq.s;
    fragment_seq.append(kseq_2->seq.s);

    string revComplement_fragment = reverse_complement(fragment_seq);


    uint64_t fragment_hash = fragmentHasher(fragment_seq);
    uint64_t revComplement_fragment_hash = fragmentHasher(revComplement_fragment);


    if (fragment_hash < revComplement_fragment_hash) {
        return fragment_hash;
    } else {
        return revComplement_fragment_hash;
    }

}


int main(int argc, char **argv) {

    if (argc != 4) {
        cerr << "run ./kDeduper <R1> <R2> <output_prefix>" << endl;
        exit(1);
    }

//    string R1_file = "/home/mabuelanin/Desktop/dev-plan/kDeduper/sample_data/R1.fastq"; // argv[1];
//    string R2_file = "/home/mabuelanin/Desktop/dev-plan/kDeduper/sample_data/R2.fastq"; // argv[2];

    string R1_file = argv[1];
    string R2_file = argv[2];
    string output_prefix = argv[3];

    if (!file_exists(R1_file)) {
        throw std::runtime_error("Could not open R1 file");
    }

    if (!file_exists(R2_file)) {
        throw std::runtime_error("Could not open R2 file");
    }

    int kSize = 31;
//    cerr << "counting number of reads ..." << endl;
//    int count = 0;
//    string line;
//    ifstream file(R1_file);
//    while (getline(file, line)) count++;
//    int no_of_sequences = count / 2;

    // --------------------------
    //       FIRST ROUND
    // ---------------------------

    gzFile fp_1, fp_2;
    kseq_t *kseq_1, *kseq_2;
    fp_1 = gzopen(R1_file.c_str(), "r");
    fp_2 = gzopen(R2_file.c_str(), "r");

    kseq_1 = kseq_init(fp_1);
    kseq_2 = kseq_init(fp_2);

    flat_hash_map<uint64_t, uint32_t> kmers_to_count;

    int terminal_kmer_offset = 2;

    for (int seqCounter = 0; kseq_read(kseq_1) >= 0 && kseq_read(kseq_2) >= 0; seqCounter++) {

        uint32_t seq_2_length = string(kseq_2->seq.s).size();

        // Extract the start and end kmers
        string start_kmer = string(kseq_1->seq.s).substr(0 + terminal_kmer_offset, kSize);
        string end_kmer = string(kseq_2->seq.s).substr(seq_2_length - kSize - terminal_kmer_offset, kSize);

        // Get the representative kmer and increment its count
        uint64_t rep_kmer = get_representative_kmer(start_kmer, end_kmer);
        kmers_to_count[rep_kmer]++;

    }

    kseq_destroy(kseq_1);
    kseq_destroy(kseq_2);
    gzclose(fp_1);
    gzclose(fp_2);


    // --------------------------
    //       SECOND ROUND
    // ---------------------------

    flat_hash_map<uint64_t, vector<uint64_t>> dedup;

    fp_1 = gzopen(R1_file.c_str(), "r");
    fp_2 = gzopen(R2_file.c_str(), "r");

    kseq_1 = kseq_init(fp_1);
    kseq_2 = kseq_init(fp_2);

    hash<string> fragmentHasher;

    fastqWriter rawReadsWriter(output_prefix + "_raw");
    fastqWriter dedupReadsWriter(output_prefix + "_dedup");


    for (int seqCounter = 0; kseq_read(kseq_1) >= 0 && kseq_read(kseq_2) >= 0; seqCounter++) {

        uint32_t seq_2_length = string(kseq_2->seq.s).size();

        // Extract the start and end kmers
        string start_kmer = string(kseq_1->seq.s).substr(0 + terminal_kmer_offset, kSize);
        string end_kmer = string(kseq_2->seq.s).substr(seq_2_length - kSize - terminal_kmer_offset, kSize);

        // Get the representative kmer and increment its count
        uint64_t rep_kmer = get_representative_kmer(start_kmer, end_kmer);

        // Get the canonical hash of the whole fragment
        uint64_t fragment_canonical_hash = canonicalFragmentHash(kseq_1, kseq_2, fragmentHasher);

        if (kmers_to_count[rep_kmer] > 1) {
            if (dedup.find(rep_kmer) == dedup.end()) {
                dedup[rep_kmer].emplace_back(fragment_canonical_hash);
                dedupReadsWriter.write(kseq_1, kseq_2);
            } else {
                if (find(dedup[rep_kmer].begin(), dedup[rep_kmer].end(), fragment_canonical_hash) ==
                    dedup[rep_kmer].end()) {
                    dedup[rep_kmer].emplace_back(fragment_canonical_hash);
                    dedupReadsWriter.write(kseq_1, kseq_2);
                }
            }
        } else {
            rawReadsWriter.write(kseq_1, kseq_2);
            dedup[rep_kmer].emplace_back(fragment_canonical_hash);
            continue;
        }

    }

}
