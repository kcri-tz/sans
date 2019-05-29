struct SANS_opt : CCDBG_Build_opt {

    string filename_colors_in;

    size_t top_splits;
    
    bool allow_asym;

    SANS_opt() :  top_splits(0), allow_asym(false) {}
};

