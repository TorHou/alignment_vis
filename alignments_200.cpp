#include <fstream>
#include <seqan/align.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <seqan/basic.h>


using namespace seqan;

int main(int argc, char const ** argv)
{
    if (argc != 3){
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
    //std::cout << id1 << std::endl;
    //std::cout << seq1 << std::endl;
    //std::cout << id2 << std::endl;
    //std::cout << seq2 << std::endl;

    int max_overlap = 200;
    

    TStringSet sequences;
    TAlign alignment;
    int score;
    
    clear(sequences);
    int len1 = length(seq1);
    int len2 = length(seq2);
    //int len1 = 5;
    resize(sequences, len1);
    
    Suffix<DnaString >::Type suf;
    int seq1start = 0;
    if (len1 >= max_overlap){
        seq1start = len1 - max_overlap;
    }
    int seq2prefixlen = len2;
    if (len2 >= max_overlap){
        seq2prefixlen = max_overlap;
    }
    Prefix<DnaString>::Type pref2 = prefix(seq2, seq2prefixlen);

    for (int i = seq1start; i < len1; ++i){
        suf = suffix(seq1, i);
        clear(sequences);
        appendValue(sequences, pref2);
        appendValue(sequences, suf);
        //std::cout << i << ": " << suf << std::endl;
        alignment = TAlign(sequences);
        score = globalAlignment(alignment, Score<int, Simple>(1, -1, -1), AlignConfig<false, false, false, true>(), LinearGaps());
        float fscore = (float) score/ (float)(len1-i);
        std::cout << fscore << std::endl;

    }

    return 0;
}
/*

        //std::cout << all_bridges[i] << std::endl;
        appendValue(sequences, seq1);
        appendValue(sequences, all_bridges[i]);
        clear(alignG);
        alignG = TAlignGraph(sequences);
        score = globalAlignment(alignG, Score<int, Simple>(2, -2, -2), AlignConfig<false, false, true, true>(), LinearGaps());
        //std::cout << "Score: " << score << std::endl;
        //std::cout << alignG << std::endl;
        int tpos = (int) getLastCoveredPosition(alignG,0);
        int nid1 = -1;
        int npos1 = -1;
        getProjectedPosition(alignG,0, tpos-1, nid1, npos1);
        //std::cout << "last covered: " << (int) 0 << "\t" << (int) tpos << std::endl;
        //std::cout << "last covered: " << (int) nid << "\t" << (int) npos << std::endl;
        Suffix<DnaString >::Type suf = suffix(all_bridges[i], npos1+1);
        //std::cout << "Suffix: " << suf << std::endl;

        //std::cout << all_bridges[i] << std::endl;
        clear(sequences);
        appendValue(sequences, seq2);
        appendValue(sequences, all_bridges[i]);
        clear(alignG);
        alignG = TAlignGraph(sequences);
        score = globalAlignment(alignG, Score<int, Simple>(2, -2, -2), AlignConfig<false, true, false, true>(), LinearGaps());
        //std::cout << "Score: " << score << std::endl;
        //std::cout << alignG << std::endl;
        tpos = (int) getFirstCoveredPosition(alignG,0);
        int nid2 = -1;
        int npos2 = -1;
        getProjectedPosition(alignG,0, tpos, nid2, npos2);
        //std::cout << "fist covered: " << (int) 0 << "\t" << (int) tpos << std::endl;
        //std::cout << "first covered: " << (int) nid2 << "\t" << (int) npos2 << std::endl;
        Prefix<DnaString >::Type pref = prefix(all_bridges[i], npos2+1);
        //std::cout << "Prefix: " << pref << std::endl;
        Infix<DnaString >::Type inf = infix(all_bridges[i], npos1+1,npos2+1);
        std::cout << "Infix: " << inf << std::endl;
        //appendValue(infixset, inf);
        assignSource(row(infixset, i), inf);

    }

    globalMsaAlignment(infixset, SimpleScore(2,-2,-2,-2));
    std::cout << "Infixe: " << infixset << std::endl;

    String<ProfileChar<Dna> > profile;
    resize(profile, length(row(infixset, 0)));
    for (unsigned rowNo = 0; rowNo < 4u; ++rowNo)
        for (unsigned i = 0; i < length(row(infixset, rowNo)); ++i)
            profile[i].count[ordValue(getValue(row(infixset, rowNo), i))] += 1;

   // call consensus from this string
    DnaString consensus;
    for (unsigned i = 0; i < length(profile); ++i)
    {
        int idx = getMaxIndex(profile[i]);
        if (idx < 4)  // is not gap
            appendValue(consensus, Dna(getMaxIndex(profile[i])));
    }

    std::cout << "consensus sequence is\n"
              << consensus << "\n";
*/

    /* this could be a good check
    int score = globalAlignment(alignG1, Score<int, Simple>(1, -1, -1), AlignConfig<true, true, true, true>(), LinearGaps());
    std::cout << "Score: " << score << std::endl;
    std::cout << alignG1 << std::endl;
    */

    /* this doesnt work
    score = globalAlignment(alignG2, Score<int, Simple>(1, -1, -1), AlignConfig<false, true, false, true>(), LinearGaps());
    std::cout << "Score: " << score << std::endl;
    std::cout << alignG2 << std::endl;
    */

    /* this works
    score = globalAlignment(alignG3, Score<int, Simple>(1, -1, -1), AlignConfig<true, false, true, false>(), LinearGaps());
    std::cout << "Score: " << score << std::endl;
    std::cout << alignG3 << std::endl;
    */


