#ifndef COMMON_H_
#define COMMON_H_

#include <seqan3/io/all.hpp>
#include <seqan3/search/all.hpp>
#include <utility>      // std::pair, std::make_pair

using namespace seqan3;


template <typename T>
struct Tag {};

struct HammingDistance_;
struct LevenshteinDistance_;

typedef Tag<HammingDistance_>       HammingDistance;
typedef Tag<LevenshteinDistance_>   EditDistance;

typedef sdsl::rank_support_v<> TSupport;

/*
template <typename TChar, typename TAlloc, typename TOwner>
struct SAValue<StringSet<String<TChar, TAlloc>, TOwner> >
{
    typedef Pair<uint16_t, uint32_t> Type;
};

template <typename TChar, typename TAlloc>
struct SAValue<String<TChar, TAlloc> >
{
    typedef uint32_t Type;
};*/

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

struct hit{
    bool rev = false;
    std::pair <uint16_t, uint32_t> occ {0,0};
    uint8_t errors = 0;
//     dna4_vector read;
};


bool occ_smaller(const hit & x, const hit & y)
{
    if(x.occ.first == y.occ.first){
        if(x.occ.second == y.occ.second)
            return x.errors < y.errors;
        else
            return x.occ.second < y.occ.second;
    }
    else
    {
        return x.occ.first < y.occ.first;
    }
}

bool occ_same(const hit & x, const hit & y)
{
    return(x.occ.first == y.occ.first && x.occ.second == y.occ.second && x.errors == y.errors);
}

template<int disT>
bool occ_similar(const hit & x, const hit & y)
{
    return(x.occ.first == y.occ.first && x.occ.second + disT >= y.occ.second && x.occ.second <= y.occ.second + disT);
}

std::string mytime()
{
    auto r = time(nullptr);
    auto c = ctime(&r);
    std::string buf(c);
    buf.insert(0, "[");
    buf.append("] ");
    buf.erase(remove(buf.begin(), buf.end(), '\n'), buf.end());
    return buf;
}

struct SearchParams
{
    unsigned length;
    unsigned overlap;
    unsigned threads;
    // bool indels;
    static constexpr bool outputProgress = false;
};

inline bool open(bi_fm_index<dna4_vector/*, my_fm_index_traits*/> & index, dna4_vector & text, const std::string & indexPath)
{
    index.load(indexPath);
    seq_data genome;
    sequence_file_input<my_traits> finindex{indexPath};
    genome.sequences = std::move(get<field::SEQ>(finindex));
    text = genome.sequences[0];
    index.text = & text;
}

inline bool file_exists(const std::string& name) {
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0);
}

//TODO fix this for multiple sequences
template <typename TIndex>
auto getSeqLengths(bi_fm_index<dna4_vector/*, my_fm_index_traits*/> & index){
//     auto mylimits = stringSetLimits(indexText(index)); //TODO test this
//     for(int i = 2; i < length(mylimits); ++i)
//         mylimits[i] -= mylimits[i - 1];
//     return(mylimits);

    std::vector<uint32_t> sl;
    sl.push_back(0);
    sl.push_back(index.fwd_fm.index.size());
    return(sl);

//     auto const & genome = indexText(index);
//     std::vector<uint32_t> sl;
//     sl.push_back(0);
//     for(uint32_t i = 0; i < countSequences(index)/*seqan::length(genome)*/; ++i)
//         sl.push_back(seqan::length(genome[i]));
//     return sl;
}


#endif