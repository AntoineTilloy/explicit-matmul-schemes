// Complete rational tensor replay, using one signed128-bit U-coordinate slice.
#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using I = __int128_t;
using U = __uint128_t;

static std::string decimal(I value) {
    bool negative = value < 0;
    U x = negative ? U(-value) : U(value);
    std::string out;
    do { out.push_back(char('0' + x % 10)); x /= 10; } while (x);
    if (negative) out.push_back('-');
    std::reverse(out.begin(), out.end());
    return out;
}

static void require(bool yes, const char *message) {
    if (!yes) throw std::runtime_error(message);
}

template<class T> static T scalar(std::ifstream &in) {
    T value{};
    in.read(reinterpret_cast<char *>(&value), sizeof(value));
    require(bool(in), "truncated input");
    return value;
}

static std::vector<int64_t> vector64(std::ifstream &in, size_t size) {
    std::vector<int64_t> value(size);
    in.read(reinterpret_cast<char *>(value.data()), size * sizeof(int64_t));
    require(bool(in), "truncated vector");
    return value;
}

struct Factor { std::vector<int64_t> ptr, idx, num; int64_t maximum = 0; };

int main(int argc, char **argv) {
    try {
        require(argc == 3, "usage: checker binary output_directory");
        std::ifstream in(argv[1], std::ios::binary);
        require(bool(in), "cannot open binary");
        std::array<char, 8> magic{};
        in.read(magic.data(), magic.size());
        require(std::string(magic.data(), magic.size()) == "ODD12801", "bad magic");
        const int64_t n = scalar<int64_t>(in), rank = scalar<int64_t>(in);
        require(n == 31 && rank == 13433, "wrong primitive dimensions/rank");
        std::array<int64_t, 3> scales;
        I scale_product = 1;
        for (auto &s : scales) { s = scalar<int64_t>(in); require(s > 0, "bad scale"); scale_product *= s; }
        require(scales == std::array<int64_t, 3>{366912, 312, 312}, "wrong pinned scales");
        const int64_t size = n * n;
        std::array<Factor, 3> factors;
        for (auto &f : factors) {
            const int64_t nnz = scalar<int64_t>(in);
            require(nnz >= 0 && nnz < 10000000, "invalid nnz");
            f.ptr = vector64(in, rank + 1);
            f.idx = vector64(in, nnz);
            f.num = vector64(in, nnz);
            require(f.ptr.front() == 0 && f.ptr.back() == nnz, "bad CSR endpoints");
            std::vector<int64_t> seen(size, -1);
            for (int64_t t = 0; t < rank; ++t) {
                require(f.ptr[t] >= 0 && f.ptr[t] <= f.ptr[t+1] && f.ptr[t+1] <= nnz, "bad CSR pointer");
                for (int64_t e = f.ptr[t]; e < f.ptr[t+1]; ++e) {
                    require(f.idx[e] >= 0 && f.idx[e] < size, "index outside factor");
                    require(seen[f.idx[e]] != t, "duplicate CSR coordinate in term");
                    seen[f.idx[e]] = t;
                    require(f.num[e] != std::numeric_limits<int64_t>::min(), "unsafe absolute value");
                    f.maximum = std::max(f.maximum, f.num[e] < 0 ? -f.num[e] : f.num[e]);
                }
            }
        }
        require(in.peek() == std::ifstream::traits_type::eof(), "unexpected trailing binary data");
        U bound = U(rank);
        const U signed_limit = (U(1) << 127) - 1;
        for (const auto &f : factors) {
            require(f.maximum > 0 && bound <= signed_limit / U(f.maximum), "overflow in bound");
            bound *= U(f.maximum);
        }
        require(U(scale_product) <= signed_limit - bound, "unsafe128-bit accumulation bound");
        bound += U(scale_product);
        const I pair_bound = I(factors[0].maximum) * factors[1].maximum;
        require(pair_bound <= std::numeric_limits<int64_t>::max(), "unsafe64-bit U/V product");
        const auto &a=factors[0], &b=factors[1], &c=factors[2];
        uint64_t expected_updates = 0;
        for (int64_t t = 0; t < rank; ++t)
            expected_updates += uint64_t(a.ptr[t+1]-a.ptr[t]) * uint64_t(b.ptr[t+1]-b.ptr[t]) * uint64_t(c.ptr[t+1]-c.ptr[t]);
        require(expected_updates == 970421136ULL, "wrong pinned expansion count");
        std::vector<size_t> offsets(size+1);
        for (int64_t x : a.idx) ++offsets[x+1];
        for (int64_t x = 0; x < size; ++x) offsets[x+1] += offsets[x];
        auto cursors = offsets;
        std::vector<int64_t> terms(a.idx.size()), nums(a.idx.size());
        for (int64_t t = 0; t < rank; ++t)
            for (int64_t e = a.ptr[t]; e < a.ptr[t+1]; ++e) {
                const size_t at = cursors[a.idx[e]]++;
                terms[at]=t; nums[at]=a.num[e];
            }
        std::vector<I> slice(size * size);
        uint64_t updates=0, scanned=0;
        const std::string progress = std::string(argv[2]) + "/odd_native_progress.json";
        for (int64_t ui = 0; ui < size; ++ui) {
            std::fill(slice.begin(), slice.end(), I(0));
            for (size_t e = offsets[ui]; e < offsets[ui+1]; ++e) {
                const int64_t t=terms[e], u=nums[e];
                for (int64_t ve=b.ptr[t]; ve<b.ptr[t+1]; ++ve) {
                    const int64_t uv=u*b.num[ve];
                    I *row=slice.data()+b.idx[ve]*size;
                    for (int64_t we=c.ptr[t]; we<c.ptr[t+1]; ++we)
                        row[c.idx[we]] += I(uv)*c.num[we];
                    updates += c.ptr[t+1]-c.ptr[t];
                }
            }
            const int64_t i=ui/n, k=ui%n;
            for (int64_t j=0; j<n; ++j) slice[(k*n+j)*size+i*n+j] -= scale_product;
            for (size_t pos=0; pos<slice.size(); ++pos) {
                ++scanned;
                if (slice[pos]) {
                    std::cout << "{\"status\":\"EXACT_TENSOR_MISMATCH\",\"u\":" << ui
                              << ",\"v\":" << pos/size << ",\"w\":" << pos%size
                              << ",\"scaled_difference\":\"" << decimal(slice[pos]) << "\"}\n";
                    return 2;
                }
            }
            if (ui % 100 == 0 || ui+1 == size) {
                std::ofstream p(progress+".tmp");
                p << "{\"completed_slices\":" << ui+1 << ",\"total_slices\":" << size
                  << ",\"updates\":" << updates << ",\"scanned\":" << scanned << "}\n";
                p.close(); require(std::rename((progress+".tmp").c_str(),progress.c_str())==0,"progress rename failed");
            }
        }
        require(updates==expected_updates && scanned==uint64_t(size)*size*size,"incomplete tensor traversal");
        std::cout << "{\"status\":\"EXACT_ODD31_FULL_TENSOR_PASS\",\"exact\":true,\"rank\":" << rank
                  << ",\"tensor_coefficients\":" << scanned << ",\"cartesian_expansions\":" << updates
                  << ",\"integer_accumulation_bound\":\"" << decimal(I(bound))
                  << "\",\"slice_bytes\":" << slice.size()*sizeof(I)
                  << ",\"CSR_unique\":true,\"native_128bit_bound_checked\":true}\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "checker error: " << e.what() << '\n';
        return 1;
    }
}
