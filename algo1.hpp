using namespace seqan3;

template <unsigned errors, typename TIndex, typename TText, typename TContainer, typename TDistanceTag>
inline void runAlgoTrivial(TIndex & index, TText const & text, TContainer & c, SearchParams const & searchParams, TDistanceTag const &)
{
/* //add errors and !needle! to delegate
    typedef typename TContainer::value_type value_type;

    constexpr uint64_t max_val = std::numeric_limits<typename TContainer::value_type>::max();
    auto scheme = OptimalSearchSchemes<0, errors>::VALUE; //THIS
    _optimalSearchSchemeComputeFixedBlocklength(scheme, searchParams.length);//THIS
    _optimalSearchSchemeComputeChronBlocklength(scheme);//THIS

    uint64_t textLength = text.size();

    #pragma omp parallel for schedule(dynamic, 1000000) num_threads(searchParams.threads)
    for (uint64_t i = 0; i < textLength - searchParams.length + 1; ++i)
    {
        value_type hits = 0;

        std::vector<hit> myhits;

//TODO add different procedures for hamming and edit distance
//         auto delegate = [&hits](auto const &it, auto const & , unsigned const ) {
//             if ((uint64_t) hits + countOccurrences(it) <= max_val)
//                 hits += countOccurrences(it);
//             else
//                 hits = max_val;
//         };
//
//
//         auto delegateDirect = [&hits](auto const & pos, DnaString const & , uint8_t const )
//         {
//             if((uint64_t) hits + 1 <= max_val)
//                 ++hits;
//         };

        auto delegate = [&myhits](auto const &it, auto const & needle, unsigned const occErrors, bool const rev) {
           for (auto occ : it.locate()){
                hit me;
                me.occ.second = occ;
//                 me.errors = occErrors;
                myhits.push_back(me);
            }
        };

        auto delegateDirect = [&myhits](auto const & pos, DnaString const & needle, uint8_t const occErrors)
        {
            hit me;
            me.occ.second = pos;
//             me.errors = occErrors;
            myhits.push_back(me);
        };


        DnaString const & needle = infix(text, i, i + searchParams.length);
        Iter<TIndex, VSTree<TopDown<> > > it(index);

        _optimalSearchScheme(delegate, delegateDirect, it, needle, scheme, EditDistance());

        std::sort(myhits.begin(), myhits.end(), occ_smaller);


//         if(i < 1){
//             for(int j = 0; j < myhits.size(); ++j)
//                 std::cout << myhits[j].occ << "  "  << (int)myhits[j].errors << "\n";
//         }

        myhits.erase(std::unique(myhits.begin(), myhits.end(), occ_similar<10>), myhits.end());

//         if(i < 1){
//             for(int j = 0; j < myhits.size(); ++j)
//                 std::cout << myhits[j].occ << "  "  << (int)myhits[j].errors << "\n";
//         }


        c[i] = myhits.size();

//         c[i] = hits;

//         std::cout << "Needle: " <<  needle << "\n";
//         std::cout << (int)hits << "\n";


    }

    resetLimits(indexText(index), c, searchParams.length);
    */

}

