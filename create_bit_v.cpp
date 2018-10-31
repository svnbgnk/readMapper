#include <iostream>
#include <seqan3/argument_parser/all.hpp>
#include <sdsl/bit_vectors.hpp>

#include "common.h"

using namespace seqan3;
using namespace std;

struct bitvectors
{
    vector<bool> fwdd;
    vector<string> names;
    vector<sdsl::bit_vector> bv;
};

vector<uint16_t> read(const string mappability_path){
    string mappability_str;

    vector<uint16_t> mappability_int;
    ifstream file(mappability_path, std::ios::binary);
    if (!file.eof() && !file.fail())
    {
        file.seekg(0, std::ios_base::end);
        std::streampos fileSize = file.tellg();
        mappability_int.resize(fileSize);
        file.seekg(0, std::ios_base::beg);
        file.read((char*)&mappability_int[0], fileSize);
        file.close();
        cout << "Load successful" << endl;
        }
    return(mappability_int);
}


template <unsigned errors>
bitvectors create_all_bit_vectors(const vector <uint16_t> & mappability, uint32_t const len, int16_t const offset, uint32_t const threshold){
    bitvectors b;
    int e = errors;
    /*
    auto scheme = OptimalSearchSchemes<0, errors>::VALUE;
    _optimalSearchSchemeComputeFixedBlocklength(scheme, len);
    _optimalSearchSchemeComputeChronBlocklength(scheme);
    auto s = scheme[0];
    sdsl::bit_vector righti (mappability.size() + len - 1, 0);
    sdsl::bit_vector lefti (mappability.size() + len - 1, 0);
    #pragma omp parallel for schedule(static)
    for(uint32_t i = 0; i < mappability.size() - offset; ++i){
        lefti[i + len - 1 + offset] = (mappability[i] <= threshold);
        righti[i + offset] = (mappability[i] <= threshold);
    }
    cout << "Finished Default Bit Vectors.  Length: " << righti.size() << endl;


    b.bv.push_back(righti);
    b.names.push_back("r_bit_vector_" + to_string(len) + "_" + to_string(e) + "_shift_0");
    b.fwdd.push_back(true);

    uint8_t blocks = s.pi.size();
    if(errors != 0){
        for(uint32_t i = 0; i < blocks - 1; ++i){
            sdsl::bit_vector newright(mappability.size() + len - 1, 0); //TODO think 0 or 1 in edge cases
            uint32_t shift = s.chronBL[i];
            cout << "r bitvector  name: " << to_string(i + 1) << endl;
            cout << "with shift: " << shift << endl;

            #pragma omp parallel for schedule(static)
            for(uint32_t j = 0; j < righti.size(); ++j){
                if(j >= shift)
                    newright[j] = righti[j - shift];
            }
            b.bv.push_back(newright);
            b.names.push_back("r_bit_vector_" + to_string(len) + "_" + to_string(e) + "_shift_" + to_string(i + 1));
            b.fwdd.push_back(true);
        }

        for(uint32_t i = 1; i < blocks; ++i){
            sdsl::bit_vector newleft(mappability.size() + len - 1, 0);//TODO think 0 or 1 in edge cases
            uint32_t shift = s.revChronBL[blocks - i];
            cout << "l bitvector  name: " << to_string(i) << endl;
            cout << "with shift: " << shift << endl;
            #pragma omp parallel for schedule(static)
            for(uint32_t j = 0; j < righti.size(); ++j){
                if(j + shift < lefti.size() - 1)
                    newleft[j] = lefti[j + shift];
            }
            b.bv.push_back(newleft);
            b.names.push_back("l_bit_vector_" + to_string(len) + "_" + to_string(e) + "_shift_" + to_string(i));
            b.fwdd.push_back(false);
        }
    }

    b.bv.push_back(lefti);
    b.names.push_back("l_bit_vector_" + to_string(len) + "_" + to_string(e) + "_shift_0");
    b.fwdd.push_back(false);
    */
    return(b);
}


template <unsigned errors>
bitvectors create_bit_vectors(const vector <uint16_t> & mappability, uint32_t const len, int16_t const offset, uint32_t const threshold){

    int e = errors;
    cout << "Create minimum amount of bitvectors" << endl;
    bitvectors b;
    sdsl::bit_vector righti (mappability.size() + len - 1, 0);
    sdsl::bit_vector lefti (mappability.size() + len - 1, 0);
    #pragma omp parallel for schedule(static)
    for(uint32_t i = 0; i < mappability.size() - offset; ++i){
        lefti[i + len - 1 + offset] = (mappability[i] <= threshold);
        righti[i + offset] = (mappability[i] <= threshold);
    }
    cout << "Finished Default Bit Vectors.  Length: " << righti.size() << endl;
    vector<sdsl::bit_vector> bit_vectors;
    vector<string> names;


    /*
    auto scheme = OptimalSearchSchemes<0, errors>::VALUE;
    _optimalSearchSchemeComputeFixedBlocklength(scheme, len);
    _optimalSearchSchemeComputeChronBlocklength(scheme);
    for (auto & s : scheme){
        bool fwd = (s.pi[0] < s.pi[1]);
        uint32_t pos = s.pi[0];
        cout << pos << endl;
        cout << "Direction forward " << fwd << endl;
        uint8_t blocks = s.pi.size();
        if(pos == 1)
            {
                b.bv.push_back(righti);
                b.names.push_back("right_bit_vector_" + to_string(len) + "_" + to_string(e));
                b.fwdd.push_back(true);
                cout << "case1" << endl;
            }
        else if(pos == blocks)
            {
                b.bv.push_back(lefti);
                b.names.push_back("left_bit_vector_" + to_string(len) + "_" + to_string(e));
                b.fwdd.push_back(false);
                cout << "case2" << endl;
            }
        else
            {
                if(fwd){
                    sdsl::bit_vector newright(mappability.size() + len - 1, 0);
                    uint32_t shift = s.chronBL[pos - 2];
                    cout << "shift r_bit for: " << shift << endl;
                    cout << "pos:  " << pos << endl;
                    for(uint32_t j = 0; j < righti.size(); ++j){
                        if(j >= shift)
                            newright[j] = righti[j - shift];
                    }
                    b.bv.push_back(newright);
                    b.names.push_back("middle_bit_vector_" + to_string(len) + "_" + to_string(e));
                    b.fwdd.push_back(true);
                }else{
                    //NOTE int32_t is not large enought to handle all the positions in the hg
                    sdsl::bit_vector newleft(mappability.size() + len - 1, 0);
                    uint32_t shift = s.revChronBL[pos];
                    cout << "shift l_bit for: " << shift << endl;
                    cout << "pos:  " << pos << endl;
                    for(uint32_t j = 0; j < righti.size(); ++j){
                        if(j + shift < lefti.size() - 1)
                            newleft[j] = lefti[j + shift];
                    }
                    b.bv.push_back(newleft);
                    b.names.push_back("middle_bit_vector_" + to_string(len) + "_" + to_string(e));
                    b.fwdd.push_back(false);
                }
            }
    }
    */
    return(b);
}

//TODO create template for mappability vector
bitvectors create_bit_vectors(const vector <uint16_t> & mappability, uint32_t const len, uint32_t const threshold, int16_t const offset, bool const bit3, uint8_t const errors){
    bitvectors result;
    if(bit3){
        switch (errors)
        {
            case 0: result = create_bit_vectors<0>(mappability, len, offset, threshold);
                    break;
            case 1: result = create_bit_vectors<1>(mappability, len, offset, threshold);
                    break;
            case 2: result = create_bit_vectors<2>(mappability, len, offset, threshold);
                    break;
            case 3: result = create_bit_vectors<3>(mappability, len, offset, threshold);
                    break;
            default: cerr << "E = " << errors << " not yet supported.\n";
                    exit(1);
        }
    }else{
        switch (errors)
        {
            case 0: result = create_all_bit_vectors<0>(mappability, len, offset, threshold);
                    break;
            case 1: result = create_all_bit_vectors<1>(mappability, len, offset, threshold);
                    break;
            case 2: result = create_all_bit_vectors<2>(mappability, len, offset, threshold);
                    break;
            case 3: result = create_all_bit_vectors<3>(mappability, len, offset, threshold);
                    break;
            default: cerr << "E = " << errors << " not yet supported.\n";
                    exit(1);
        }
    }
    return(result);
}


void order_bit_vector(bitvectors & b, std::string const indexPath, uint32_t const threads)
{
    /*
    cout << mytime() << "Loading Index" << endl;

    bi_fm_index<dna4_vector> index; //my_fm_index_traits
    dna4_vector text;
    open(index, text, indexPath);

    uint32_t indexSize = index.fwd_fm.index.size(); //TODO fix this
    cout << mytime() << "Loaded Index. Size:" << indexSize << endl;
    vector<sdsl::bit_vector> bit_vectors_ordered (b.bv);
    uint32_t number_of_indeces = 1;//countSequences(index); //TODO fix for multiple sequences

    std::vector<uint32_t> sequenceLengths = getSeqLengths(index);
    for(uint32_t i = 2; i < sequenceLengths.size(); ++i)
        sequenceLengths[i] += sequenceLengths[i - 1];

    cout << "Number of Sequences in index: " << 1 << endl; //countSequences(index)
    cout << mytime() << "Start sorting bitvectors" << endl;

    uint32_t mythreads;
    if(threads == 0)
        mythreads = 1 ; //omp_get_max_threads()
    else
        mythreads = threads;
    //dynamic since tasks can take different amount of time (SA Sampling) critical to guarantee thread safety

    uint8_t bsize = b.bv.size();

//     #pragma omp parallel for schedule(dynamic) num_threads(mythreads)
    #pragma omp parallel for schedule(static, (indexSize/(mythreads*100))) num_threads(mythreads)
    for (unsigned j = 0; j < indexSize - number_of_indeces; ++j)
    {
        // skip sentinels
        std::pair<uint16_t, uint32_t> sa_f = index.fwd.sa[j + number_of_indeces];
        std::pair<uint16_t, uint32_t> sa_r = index.rev.sa[j + number_of_indeces];

        uint32_t fpos = sa_f.i2 + sequenceLengths[sa_f.i1];
        uint32_t rpos = sequenceLengths[sa_r.i1 + 1] - sa_r.i2 - 1;
        vector<bool> values(bsize);

        for(uint32_t i = 0; i < bsize; ++i){
            if(b.fwdd[i]){
                values[i] = b.bv[i][fpos];
            }
            else
            {
                values[i] = b.bv[i][rpos];
            }
        }
        #pragma omp critical
        {
        for(uint32_t i = 0; i < bsize; ++i)
            bit_vectors_ordered[i][j] = values[i];
        }
    }
    b.bv = bit_vectors_ordered;
    */
}

int main(int argc, char *argv[])
{

    std::string indexPath, mappabilityPath;
    int length, errors = 0, r = 0, threshold = 10, offset = 0, threads = 1; //benchparams
    bool debug, bit3version; //stats

    argument_parser parser("Search", argc, argv);
    parser.add_option(indexPath, 'I', "index", "Path to the index directory + name"/*, , detail::default_validator*/);
    parser.add_option(mappabilityPath, 'M', "mappability", "Path to mappability file");


    parser.add_option(length, 'K', "length", "Length of reads");
    parser.add_option(errors, 'E', "errors", "Max errors allowed during mapping");
    parser.add_option(threshold, 'T', "threshold", "Number of times a k-mer can occure (needed for compare)");
    parser.add_option(offset, 'o', "offset", "Shift bitvectors to the right with positive number (needed for Edit Distance)");

    parser.add_flag(debug, 'd', "debug", "Also create chronical bit_vectors (for debugging)");
    parser.add_flag(bit3version, 'm', "bit3version", "Only create the required 3 bitvectors needed for acquiring mappability in non-unidirectional cases");
    parser.add_option(threads, 't', "threads", "Number of threads used");

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

    cout << mytime() << "Program start." << endl;
    vector<uint16_t> mappability = read(mappabilityPath);
    cout << mytime() << "Loaded Mappability vector. Size: " << mappability.size() << endl;
/*
    bitvectors result = create_bit_vectors(mappability, len, offset, threshold, bit3, errors);

    cout << mytime() << "Finished bit vectors." << endl;

    if(debug)
    {
        for(uint32_t i = 0; i < result.bv.size(); ++i){
            sdsl::store_to_file(result.bv[i], toCString(outputPath) + result.names[i] + "_chron");
            std::ofstream outfile((toCString(outputPath) + result.names[i] + "_debug"), std::ios::out | std::ofstream::binary);
            std::copy(result.bv[i].begin(), result.bv[i].end(), std::ostream_iterator<bool>(outfile));
            outfile.close();
        }
    }

    //order bitvectors according to the forward or reverse suffix array
    order_bit_vector(result, indexPath, mmap, alphabet, threads);
    cout << mytime() << "Finished sorting" << endl;
    for(uint32_t i = 0; i < result.bv.size(); ++i){
        sdsl::store_to_file(result.bv[i], toCString(outputPath) + result.names[i]);
        if(debug){
            std::ofstream outfile((toCString(outputPath) + result.names[i] + "_osa_debug"), std::ios::out | std::ofstream::binary);
            std::copy(result.bv[i].begin(), result.bv[i].end(), std::ostream_iterator<bool>(outfile));
            outfile.close();
        }
    }*/

    cout << mytime() << "Finished saving bit vectors" << endl;

    return 0;
}