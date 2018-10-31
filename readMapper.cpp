#include <iostream>
#include <thread>
#include <chrono>

#include <seqan3/argument_parser/all.hpp>
// #include <sdsl/suffix_trees.hpp>
#include <seqan3/core/all.hpp>
#include <seqan3/io/all.hpp>
#include <seqan3/search/all.hpp>
#include <seqan3/alphabet/all.hpp>
#include <seqan3/range/all.hpp>
#include <seqan3/std/view/reverse.hpp>
// #include <seqan3/io/stream/debug_stream.hpp>

#include "common.h"


using namespace seqan3;
using namespace std;
using namespace literal;

typedef sdsl::bit_vector TBitvector;
typedef sdsl::rank_support_v<> TSupport;
/*
struct my_traits : sequence_file_input_default_traits_dna
{
        using sequence_alphabet = dna4;                        // instead of dna5
//         template <typename alph>
//         using sequence_container = bitcompressed_vector<alph>; // must be defined as a template!
// what is this ???
};

struct seq_data
{
    concatenated_sequences<dna4_vector>  sequences;
    concatenated_sequences<std::string>  ids;
};

struct my_fm_index_traits
{
    // Type of the underlying SDSL index.
    using sdsl_index_type = sdsl::csa_wt<
        sdsl::wt_blcd<
            sdsl::bit_vector,
            sdsl::rank_support_v<>,
            sdsl::select_support_scan<>,
            sdsl::select_support_scan<0>
        >,
        10,
        10000000,
        sdsl::sa_order_sa_sampling<>,
        sdsl::isa_sampling<>,
        sdsl::plain_byte_alphabet
    >;
};

struct my_bi_fm_index_traits
{
    using fm_index_traits = my_fm_index_traits;

    using rev_fm_index_traits = my_fm_index_traits;
};*/

dna4_vector reverseComplement(dna4_vector const & seq)
{
    dna4_vector rev = view::reverse(seq);
    dna4_vector rc = view::complement(rev);
    return(rc);
}

seq_data createRCReads(seq_data const & reads)
{
    //TODO fix this copy mess
    seq_data rcReads;
    for(int i = 0; i < reads.ids.size(); ++i)
    {
        dna4_vector rc = reverseComplement(reads.sequences[i]);
        rcReads.sequences.push_back(rc);
        rcReads.ids.push_back(reads.ids[i]);
    }
    return(rcReads);
}

int main(int argc, char *argv[])
{

    std::string indexPath, readsPath, bitvectorPath, outputPath;
    int length, errors = 0, r = 0, threshold = 10; //benchparams
    bool ossMAP, oss, ossITV, edit, rc, sp, fr, ecompare; //stats

    argument_parser parser("Search", argc, argv);
    parser.add_option(indexPath, 'I', "index", "Path to the index directory + name"/*, , detail::default_validator*/);
    parser.add_option(readsPath, 'R', "reads", "Path to the read fasta");
    parser.add_option(bitvectorPath, 'B', "bitvectors", "Path to the bitvector directory");
    parser.add_option(outputPath, 'O', "output", "Path to output directory");


    parser.add_option(length, 'K', "length", "Length of reads");
    parser.add_option(errors, 'E', "errors", "Max errors allowed during mapping");

    bool tmp;
    parser.add_flag(tmp, 'n', "notmy", "Do NOT use search with mappability");
    ossMAP = !tmp;
    parser.add_flag(oss, 'd', "default", "Search with Optimum Search Schemes");
    parser.add_flag(ossITV, 't', "defaultT", "Search with Optimum Search Schemes with In Text Verification");
    parser.add_flag(edit, 'e', "editDistance ", "Search Edit Distance, correct bitvectors have to be selected");
    parser.add_flag(rc, 'b', "reverseComplement", "Search on both strands");
    parser.add_flag(sp, 's', "splitReads", "Split the reads in the middle");
    parser.add_flag(fr, 'f', "filterReads", "Create fastas with filtered reads in output directory");
    parser.add_flag(ecompare, 'c', "compare", "Compare Search Results");
    parser.add_option(threshold, 'T', "threshold", "Number of times a k-mer can occure (needed for compare)");


    try
    {
        parser.parse();
    }
    catch (seqan3::parser_invalid_argument const & ext) // the user did something wrong
    {
        std::cerr << "[PARSER ERROR] " << ext.what() << "\n"; // customize your error message
        return -1;
    }
    catch (seqan3::parser_interruption const &) // expected behaviour on special requests (e.g. `--help`)
    {
        return 0;
    }

    //load index
    cout << "Loading Index: " << indexPath << "\n";
    bi_fm_index<dna4_vector> myindex; //my_fm_index_traits
    myindex.load(indexPath);

    //load text
    std::string indexText{indexPath};
    indexText += ".fasta";
    sequence_file_input<my_traits> finIndex{"indexText"};
    seq_data text;
    text.sequences = std::move(get<field::SEQ>(finIndex));
    std::vector<dna4> genome = text.sequences[0];

    myindex.text = & genome;

    //load reads
    sequence_file_input<my_traits> finReads{readsPath};
    seq_data reads;

//  reads.reserve(1000);
//  reads.concat_reserve(1'000 * 1'000);

    reads.sequences = std::move(get<field::SEQ>(finReads)); // we move the buffer directly into our storage
    reads.ids = std::move(get<field::ID>(finReads)); // we move the buffer directly into our storage

    uint32_t nreads = reads.sequences.size();
    if(r > nreads)
        cerr << "Not enough reads" << endl;

//     myvector.erase (myvector.begin(),myvector.begin()+3);
    if(r != 0){
        reads.sequences.erase(reads.sequences.begin(), reads.sequences.begin() + r);
    }


    if(sp){
        auto startT = std::chrono::high_resolution_clock::now();
        int readLength = reads.sequences[0].size();
        for(int i = 0; i < reads.sequences.size(); ++i){
//             appendValue(rcReads, infix(reads[i], readlength/2 + readlength));
//             reads.sequences[i] = infix(reads[i], 0, 0 + readlength/2); //TODO what is infix?
        }
        cout << "Read length: ";
        cout << reads.sequences.size() << endl;
        auto finishT = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsedT = finishT - startT;
        cout << "Splitted reads in " << elapsedT.count() << "s" << endl;
    }

    if(rc){
        //add reverse Complement of reads to search on both strands
        auto startT = std::chrono::high_resolution_clock::now();
        seq_data rcreads = createRCReads(reads);
        reads.sequences.insert(reads.sequences.begin(), rcreads.sequences.begin(), rcreads.sequences.end());
        reads.ids.insert(reads.ids.begin(), rcreads.ids.begin(), rcreads.ids.end());
        auto finishT = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsedT = finishT - startT;
        cout << "added reversed Reads in " << elapsedT.count() << "s" << endl;
    }

    //load bitvectors
    vector<pair<TBitvector, TSupport>> bitvectors;
    if(ossMAP){
        cout << "Loading bitvectors" << endl;
//         bitvectors = loadBitvectors(bitvectorpath, K, nerrors); //TODO Implement this
//         cout << "Bit vectors loaded. Number: " << bitvectors.size() << " Length: " << bitvectors[0].first.size() << endl;
    }

    auto start = std::chrono::high_resolution_clock::now();
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed;

/*
    sequence_file_input<my_traits> finindex{"tmp/myindexText.fasta"};
    data_storage_t mytext;
    mytext.sequences = std::move(get<field::SEQ>(finindex));
    std::vector<dna4> genome2 = mytext.sequences[0];
    debug_stream << genome2 << endl;*/
    //load reads

    std::cout << ossMAP << std::endl;
    std::cout << readsPath << std::endl;
    std::cout << errors << std::endl;
//     std::cout << indexpath2 << std::endl;

/*
    std::vector<dna4> genome
    {"TGTTGTCATTGGCGGAAAACTTCCGTTCAGGAGGCGGACACTGATTGACACGGTTTAGCAGAAGGTTTGAGGAATAGGTTAAATTGAGTGGTTTAATAACGGTATGTCTGGGATTAAAGTGTAGTATAGTGTGATTATCGGAGACGGTTTTAAGACACGAGTTCCCAAAATCAAGCGGGGTCATTACAACGGTTATTCCTGGTAGTTTAGGTGTACAATGTCCTGAAGAATATTTAAGAAAAAAGCACCCCTCATCGCCTAGAATTACCTACTACGGTCGACCATACCTTCGATTATCGCGGCCACTCTCGCATTAGTCGGCAGAGGTGGTTGTGTTGCGATAGCCCAGTATAATATTCTAAGGCGTTACCCTGATGAATATCCAACGGAATTGCTATAGGCCTTGAACGCTACACGGACGATACGAAATTATGTATGGACCGGGTCATCAAAAGGTTATACCCTTGTAGTTAACATGTAGCCCGGCCCTATTAGTACAGTAGTGCCTTGAATGGCATTCTCTTTATTAAGTTTTCTCTACAGCTAAACGATCAAGTGCACTTCCACAGAGCGCGGTGGAGATTCATTCACTCGGCAGCTCTGTAATAGGGACTAAAAAAGTGATGATAATCATGAGTGCCGCGTTATGGTGGTGTCGGAACAGAGCGGTCTTACGGCCAGTCGTATGCCTTCTCGAGTTCCGTCCAGTTAAGCGTGACAGTCCCAGTGTACCCACAAACCGTGATGGCTGTGCTTGGAGTCAATCGCAAGTAGGATGGTCTCCAGACACCGGGGCACCAGTTTTCACGCCGAAAGCATAAACGACGAGCAGATATGAAAGTGTTAGAACTGGACGTGCCGTTTCTCTGCGAAGAACACCTCGAGCTGTAGCGTTGTTGCGCTGCCTAGATGCAGTGTTGCTCATATCACATTTGCTTCAACGACTGCCGCCTTCGCTGTATCCCTAGACACTCAACAGTAAGCGCTTTTTGTAGGCAGG"_dna4};
    std::vector<dna4> text{"ACGT"_dna4};
    debug_stream << text << "\n";
    concatenated_sequences<dna4_vector> reads{text};
    cout << "Here comes the debug_stream" << "\n";
    debug_stream << reads[0] << '\n';
    reads.push_back(text);
    reads.push_back(text);
// foobar.insert(foobar.end(), "ACGT"_dna4);


    cout << "Path: " << readsPath << endl;


//     sequence_file_input fin{readsPath};
//
//     for (auto & rec : fin)
//     {
//         std::cout << "ID:  " << get<field::ID>(rec) << '\n';
//         std::cout << "SEQ: " << (get<dna4_vector>(rec) | view::to_char) << '\n'; //seq retrieved is dna5 for some reason
//
//     }

    sequence_file_input<my_traits> fin{readsPath};

    data_storage_t data_storage;
    data_storage.sequences = std::move(get<field::SEQ>(fin)); // we move the buffer directly into our storage
    data_storage.ids = std::move(get<field::ID>(fin)); // we move the buffer directly into our storage

    cout << "Size: " << data_storage.sequences.size() << endl;
    cout << "Size0: " << data_storage.sequences[0].size() << endl;

     for(int i = 0; i < data_storage.sequences.size(); ++i){
        debug_stream << data_storage.ids[i] << endl << data_storage.sequences[i] << endl;
    }

    sequence_file_output fout{"my.fasta"};

    for(int i = 0; i < reads.size(); ++i){
        std::string id = to_string(i);
        dna4_vector seq = reads[i];

        fout.push_back(std::tie(seq, id));

    }
//     std::vector<dna4> text{"ACGTT"_dna4};
    dna4_vector rev = view::reverse(text);
    dna4_vector myrc = view::complement(rev);
    debug_stream << myrc << "\n";


    //     fout.close();
//     fout = range;
//     fout = std::tie(data_storage.sequences, data_storage.ids);


//     Index<StringSet<TString, Owner<ConcatDirect<> > >, TIndexConfig> index(chromosomesConcat);
    fm_index index{genome};
    auto it = index.begin();                                  // create an iterator
    it.extend_right("AAGG"_dna4);

    for (auto const & pos : it.locate())                      // outputs: 8, 22
        debug_stream << pos << ' ';
    debug_stream << '\n';


    bi_fm_index<dna4_vector> myindex{genome};
//     bi_fm_index.construct();
    auto iter = myindex.begin();
    iter.extend_right("AAC"_dna4);         // search the sequence "AAC"
    debug_stream << iter.query() << '\n';     // outputs "AAC"
    debug_stream << iter.last_char() << '\n'; // outputs 'C'
    iter.cycle_back();
    iter.extend_left(dna4::G);             // search the sequence "GAAT"
    debug_stream << iter.query() << '\n';     // outputs "GAAC"
    debug_stream << iter.last_char() << '\n'; // outputs 'G'
//     auto uni_it = iter.to_fwd_iterator(); //does forward search extend_right and cycle_back
//      auto uni_it = iter.to_rev_iterator(); //does backward search //can only extend_right cycle_back


    //TODO try load and save
    myindex.store("tmp/myindex");
//     sequence_file_output foutindex{"tmp/myindexText.fasta"};

//     data_storage_t mytextout;
//     mytextout.sequences.push_back(genome);

//     mytextout.ids.push_back(data_storage.ids[0]);
//     foutindex = std::tie(mytextout.sequences, mytextout.ids);

//      foutindex = std::tie(genome, "genome");


//     std::this_thread::sleep_for (std::chrono::seconds(1));


    sequence_file_input<my_traits> finindex{"tmp/myindexText.fasta"};
    data_storage_t mytext;
    mytext.sequences = std::move(get<field::SEQ>(finindex));
    std::vector<dna4> genome2 = mytext.sequences[0];
    debug_stream << genome2 << endl;







    bi_fm_index<dna4_vector> new_index;
    new_index.load("tmp/myindex");

    new_index.text = & genome2;
    iter = myindex.begin();
    iter.extend_right("AAC"_dna4);         // search the sequence "AAC"
    debug_stream << iter.query() << '\n';






    std::cout << "Hello SeqAn 2" << std::endl;

*/

    return 0;
}

// using TMyFastConfig = seqan::FastFMIndexConfig<void, uint32_t, 2, 1>;
// using TIndexConfig = seqan::BidirectionalIndex<seqan::FMIndex<void, TMyFastConfig> >;