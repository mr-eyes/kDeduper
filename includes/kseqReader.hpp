#include <zlib.h>
#include <cstdio>
#include "kseq.h"
#include <iostream>
#include <vector>

KSEQ_INIT(gzFile, gzread)

using namespace std;

struct peRead {
    string R1_seq;
    string R2_seq;
    string R1_name;
    string R2_name;
};


class kseqReader {

private:
    bool END = false;
    string R1_fastx, R2_fastx;
    int chunk_size;
    gzFile fp_r1, fp_r2;
    kseq_t *seq_1, *seq_2;

public:

    std::vector<peRead> chunk_reads;

    kseqReader(string R1File, string R2File, int chunk_size);

    std::vector<peRead> *next_chunk();

    [[nodiscard]] bool end() const {
        return this->END;
    }

    ~kseqReader() {
        kseq_destroy(this->seq_1);
        kseq_destroy(this->seq_2);

        gzclose(this->fp_r1);
        gzclose(this->fp_r2);
    }

};
