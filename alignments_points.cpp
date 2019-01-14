#include <stdio.h>
#include <fstream>
#include <seqan/align.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <seqan/basic.h>


using namespace seqan;

int main(int argc, char const ** argv)
{
    if (argc != 4){
        std::cout << "Wrong number of inputs, please provide two paths to the input fasta files" << std::endl;
        return 1;  // Invalid number of arguments.
    }

    typedef String<Dna> TSequence;                              // sequence type
    typedef StringSet<DnaString> TStringSet;                    // container for strings
    typedef StringSet<TSequence, Dependent<> > TDepStringSet;   // dependent string set
    typedef Align<TSequence, ArrayGaps> TAlign;                 // align type
    typedef StringSet<Infix<DnaString>> TInfixSet;

    TSequence seq1;
    TSequence seq2;
    CharString id1;
    CharString id2;
    SeqFileIn seqFileIn1(toCString(argv[1]));
    readRecord(id1, seq1, seqFileIn1);
    SeqFileIn seqFileIn2(toCString(argv[2]));
    readRecord(id2, seq2, seqFileIn2);
    int len1 = length(seq1);
    int len2 = length(seq2);
    
    int point = std::stoi(argv[3]);

    if(len1 + point < 0 or point > 0){
        std::cout << "Point is not on the first sequence.";
        exit(EXIT_FAILURE);
    }

    int start = 0;
    int end = 0;
    if (len1 + point - 10 < 0){
        start = 0;
    }
    else{
        start = len1 + point - 10;
    }
    
    if (len1 + point + 10 >= len1){
        end = len1 - 1;
    }
    else{
        end = len1 + point +10;
    }

    TStringSet sequences;
    TAlign alignment;
    int score;
    
    clear(sequences);
    resize(sequences, len1);
    
    Suffix<DnaString >::Type suf;
    int seq1start = 0;
    int seq2prefixlen = len2;

    Prefix<DnaString>::Type pref2 = prefix(seq2, seq2prefixlen);

    for (int i = start; i < end; ++i){
        suf = suffix(seq1, i);
        clear(sequences);
        appendValue(sequences, pref2);
        appendValue(sequences, suf);
        //std::cout << i << ": " << suf << std::endl;
        alignment = TAlign(sequences);
        score = globalAlignment(alignment, Score<int, Simple>(1, -1, -1), AlignConfig<false, false, false, true>(), LinearGaps());
        float fscore = (float) score/ (float)(len1-i);
        std::cout << i << "\t" << fscore << std::endl;
    }

    return 0;
}


