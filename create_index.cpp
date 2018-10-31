#include <seqan3/io/all.hpp>
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/search/all.hpp>

#include "common.h"

using namespace seqan3;
using namespace std;


// template <typename Tseq_data, typename inner_type, typename data_delimiters_type = std::vector<typename inner_type::size_type>>
// void construct_index(Tseq_data<innertype, data_delimiters_type> const & reference, std::string indexPath)
template< typename TSeqData>
void construct_index(TSeqData const & reference, std::string indexPath)
{
    //create Index
    dna4_vector genome = reference.sequences[0]; //TODO fix this
    bi_fm_index<dna4_vector/*, my_fm_index_traits*/> myindex{genome}; //TODO works only on one sequence
    myindex.store(indexPath);

    //store text
    std::string indexTextPath{indexPath};
    indexTextPath += ".fasta";
    sequence_file_output fout{indexTextPath};
    fout = std::tie(reference.sequences, reference.ids);
}

int main(int argc, char *argv[])
{
    std::string referencePath, indexPath;
    argument_parser parser("IndexCreation", argc, argv);
    parser.info.short_description = "App for creating an index. Only supports Dna (with and without N's). At most 4 gigabases in total allowed.";

    parser.add_option(referencePath, 'G', "genome", "Path to the genome");
    parser.add_option(indexPath, 'I', "index", "Path to the index directory + name");


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


    cout << "Finished Parsing" << endl;

    cout << "referencePath: " << referencePath << endl;
    cout << "indexPath: " << indexPath << endl;

    sequence_file_input<my_traits> fin{referencePath};
    seq_data reference;

    cout << "Start loading" << endl;
    reference.sequences = std::move(get<field::SEQ>(fin));
    reference.ids = std::move(get<field::ID>(fin));


    std::cout << "Number of sequences: " << reference.sequences.size() << '\n' << std::flush;
//     std::cout << "Sampling rate: " << my_fm_index_traits << '\n'; //TODO make it work this

    construct_index(reference, indexPath);

    std::cout << "Index created successfully.\n";



    return 0;
}