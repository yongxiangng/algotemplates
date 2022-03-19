struct special_hash {
    inline size_t operator()(const pair<int,vi> & v) const {
        auto &vec = v.second;
        size_t seed = vec.size();
        for(auto& i : vec) {
            seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed*31+v.first;
    }
};

