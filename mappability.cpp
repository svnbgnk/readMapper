#include <seqan3/io/all.hpp>
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/search/all.hpp>

#include "common.h"
#include "algo1.hpp"

using namespace std;
using namespace seqan3;

struct Options
{
    unsigned errors;
    bool mmap;
    bool indels;
    bool high;
    std::string indexPath;
    std::string outputPath;
    std::string alphabet;
};


string get_output_path(Options const & opt, SearchParams const & searchParams)
{
    string output_path = opt.outputPath;
    output_path += "_" + to_string(opt.errors) + "_" + to_string(searchParams.length) + "_" + to_string(searchParams.overlap);
//     output_path += ".gmapp" + string(opt.high ? "16" : "8");
    return output_path;
}

template <typename TDistance, typename value_type, typename TIndex, typename TText>
inline void run(TIndex & index, TText const & text, Options const & opt, SearchParams const & searchParams)
{
    vector<value_type> c(text.size() - searchParams.length + 1, 0);
    switch (opt.errors)
    {
        case 0:  runAlgoTrivial<0>(index, text, c, searchParams, TDistance());
                break;
        case 1:  runAlgoTrivial<1>(index, text, c, searchParams, TDistance());
                break;
        case 2:  runAlgoTrivial<2>(index, text, c, searchParams, TDistance());
                break;
        case 3:  runAlgoTrivial<3>(index, text, c, searchParams, TDistance());
                break;
        case 4:  runAlgoTrivial<4>(index, text, c, searchParams, TDistance());
                break;
        default: cerr << "E = " << opt.errors << " not yet supported.\n";
                exit(1);
    }

//     std::is_same<TDistanceTag, EditDistance>::value
/*


    if(opt.indels){
        switch (opt.errors)
        {
            case 0:  runAlgoTrivial<0>(index, text, c, searchParams);
                    break;
            case 1:  runAlgoTrivial<1>(index, text, c, searchParams);
                    break;
            case 2:  runAlgoTrivial<2>(index, text, c, searchParams);
                    break;
            case 3:  runAlgoTrivial<3>(index, text, c, searchParams);
                    break;
            case 4:  runAlgoTrivial<4>(index, text, c, searchParams);
                    break;
            default: cerr << "E = " << opt.errors << " not yet supported.\n";
                    exit(1);
        }
    }
    else
    {
    switch (opt.errors)
    {
        case 0:  runAlgo4<0>(index, text, c, searchParams);
                 break;
        case 1:  runAlgo4<1>(index, text, c, searchParams);
                 break;
        case 2:  runAlgo4<2>(index, text, c, searchParams);
                 break;
        case 3:  runAlgo4<3>(index, text, c, searchParams);
                 break;
        case 4:  runAlgo4<4>(index, text, c, searchParams);
                 break;
        default: cerr << "E = " << opt.errors << " not yet supported.\n";
                 exit(1);
    }
    }

    if (SearchParams::outputProgress)
        std::cout << '\r';
    std::cout << "Progress: 100.00%\n" << std::flush;
    cout.flush();


    string output_path = get_output_path(opt, searchParams);
    save(c, output_path);*/

}



template <typename TDistance, typename value_type>
inline void run(Options const & opt, SearchParams const & searchParams)
{
    //loadIndex
    //maybe create template for dna4 and dna5
    bi_fm_index<dna4_vector> index;
    dna4_vector text;
    open(index, text, opt.indexPath);

//     open(index, toCString(opt.indexPath), OPEN_RDONLY);
    run<TDistance, value_type>(index, text, opt, searchParams);
}



template <typename TDistance>
inline void run(Options const & opt, SearchParams const & searchParams)
{
    if (opt.high) {
        run<TDistance, uint16_t>(opt, searchParams);
    }
    else
        run<TDistance, uint8_t>(opt, searchParams);
}



// template <typename TChar, typename TAllocConfig>
inline void run(Options const & opt, SearchParams const & searchParams)
{
    if (opt.indels) {
        run<EditDistance>(opt, searchParams); //EditDistance
    }
    else
        run<HammingDistance>(opt, searchParams); //HammingDistance
}

int main(int argc, char *argv[])
{
    std::string indexPath, outputPath;
    int length, errors = 0, overlap = 0, threads = 1/*omp_get_max_threads()*/;
    bool indels, hi;

    argument_parser parser("Mappability", argc, argv);

    parser.info.short_description = "App for calculating the mappability values. Only supports Dna4 so far.";

    parser.add_option(indexPath, 'I', "index", "Path to the index directory + name"/*, , detail::default_validator*/);
    parser.add_option(outputPath, 'O', "output", "Path to output directory");


    parser.add_option(length, 'K', "length", "Length of reads");
    parser.add_option(errors, 'E', "errors", "Max errors allowed during mapping");
    parser.add_option(overlap, 'o', "overlap", "Number of overlapping reads (o + 1 Strings will be searched at once beginning with their overlap region)");
    parser.add_option(threads, 't', "threads", "Number of threads");

    parser.add_flag(indels, 'i', "indels", "Turns on indels (EditDistance)");
    parser.add_flag(hi, 'h', "reverseComplement", "Stores the mappability vector in 16 bit unsigned integers instead of 8 bit (max. value 65535 instead of 255)");

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

    // Retrieve input parameters
    Options opt;
    SearchParams searchParams;
    opt.errors = errors;
    opt.indexPath = indexPath;
    opt.outputPath = outputPath;
    opt.indels = indels;
    opt.high = hi;

    searchParams.length = length;
    searchParams.threads = threads;
    searchParams.overlap = overlap;

    if (searchParams.overlap > searchParams.length - 1)
    {
        cerr << "ERROR: overlap cannot be larger than K - 1.\n";
        exit(1);
    }

    if (!(searchParams.length - searchParams.overlap >= opt.errors + 2))
    {
        cerr << "ERROR: overlap should be at least K - E - 2. (K - O >= E + 2 must hold since common overlap has length K - O and will be split into E + 2 parts).\n";
        exit(1);
    }

    // searchParams.overlap - length of common overlap
    searchParams.overlap = searchParams.length - searchParams.overlap;




    run(opt, searchParams);




/*
    index.load(indexPath);
    seq_data texts;
    sequence_file_input<my_traits> finindex{indexPath};
    texts.sequences = std::move(get<field::SEQ>(finindex));
    dna4_vector genome = texts.sequences[0];
    index.text = & genome;*/


}