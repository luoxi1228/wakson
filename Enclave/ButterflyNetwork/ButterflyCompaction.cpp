#include "ButterflyCompaction.hpp"
#include <cstring>
#include <cstdio>

struct BNPacket
{
    uint32_t value;
    uint32_t tag;
};

static_assert(sizeof(BNPacket) == 8, "BNPacket must be 8 bytes");

void generateControlBits(const bool *selected_list, size_t n, std::vector<uint8_t> &ctrlBits)
{
    FOAV_SAFE_CNTXT(Butterfly_CONTROLBIT, n)

    ctrlBits.clear();
    if (!selected_list || n == 0)
        return;

    const uint32_t Npad = nextPow2_u32((uint32_t)n);
    if (Npad < 2u)
        return;

    uint32_t l = 0;
    for (uint32_t t = Npad; t > 1; t >>= 1)
        ++l;

    std::vector<BNPacket> P(Npad);
    for (uint32_t i = 0; i < Npad; ++i)
    {
        P[i].value = (i < n && selected_list[i]) ? 1u : 0u;
        P[i].tag = 0;
    }

    uint32_t rank = 0;
    for (uint32_t i = 0; i < Npad; ++i)
    {
        const uint32_t v = P[i].value;
        P[i].tag = rank * v;
        rank += v;
    }

    const size_t switchesPerStage = Npad / 2;
    ctrlBits.resize(switchesPerStage * l);

    size_t t = 0;
    for (uint32_t d = 0; d < l; ++d)
    {
        const uint32_t stride = (1u << d);
        for (uint32_t b = 0; b < Npad; b += (stride << 1))
        {
            for (uint32_t k = 0; k < stride; ++k)
            {
                const uint32_t i = b + k;
                const uint32_t j = b + k + stride;

                const uint32_t bit_i = (P[i].tag >> d) & 1u;
                const uint32_t bit_j = (P[j].tag >> d) & 1u;

                const uint8_t c =
                    (uint8_t)((P[i].value & bit_i) |
                              (P[j].value & (bit_j ^ 1u)));

                ctrlBits[t++] = c;

                oswap_buffer<OSWAP_8>(
                    (unsigned char *)&P[i],
                    (unsigned char *)&P[j],
                    8u,
                    c);
            }
        }
    }
}

static inline void fillRandomSelected(bool *selected_list, size_t N)
{
    for (size_t i = 0; i < N; ++i)
    {
        uint8_t r = 0;
        getRandomBytes(&r, 1);
        selected_list[i] = (bool)(r & 1u);
    }
}

static inline void applyButterflyCompact_dispatch(unsigned char *buf, size_t N, size_t block_size, const std::vector<uint8_t> &ctrlBits)
{
    if (block_size == 4)
    {
        applyCompaction<OSWAP_4>(buf, N, block_size, ctrlBits);
    }
    else if (block_size == 8)
    {
        applyCompaction<OSWAP_8>(buf, N, block_size, ctrlBits);
    }
    else if (block_size == 12)
    {
        applyCompaction<OSWAP_12>(buf, N, block_size, ctrlBits);
    }
    else if (block_size % 16 == 0)
    {
        applyCompaction<OSWAP_16X>(buf, N, block_size, ctrlBits);
    }
    else
    {
        applyCompaction<OSWAP_8_16X>(buf, N, block_size, ctrlBits);
    }
}

double testButterflyCompaction(unsigned char *buffer, size_t N, size_t block_size, bool *selected_list, enc_ret *ret)
{
    FOAV_SAFE2_CNTXT(testButterflyCompaction, N, block_size)

    if (!buffer || N < 2)
    {
        if (ret)
        {
            ret->ptime = 0.0;
            ret->control_bits_time = 0.0;
            ret->apply_perm_time = 0.0;
            ret->OSWAP_count = 0;
            ret->OSWAP_cb = 0;
            ret->OSWAP_ap = 0;
        }
        return 0.0;
    }

    long t1, t2;

    // ----------------
    // Offline phase
    // ----------------
#ifdef COUNT_OSWAPS
    OSWAP_COUNTER = 0;
#endif

    std::vector<uint8_t> ctrlBits;
    ocall_clock(&t1);
    generateControlBits(selected_list, N, ctrlBits);
    ocall_clock(&t2);

    double cb_ms = ((double)(t2 - t1)) / 1000.0;

#ifdef COUNT_OSWAPS
    size_t oswap_cb = OSWAP_COUNTER;
#else
    size_t oswap_cb = 0;
#endif

    // ----------------
    // Online phase
    // ----------------
#ifdef COUNT_OSWAPS
    OSWAP_COUNTER = 0;
#endif

    ocall_clock(&t1);
    applyButterflyCompact_dispatch(buffer, N, block_size, ctrlBits);
    ocall_clock(&t2);

    double ap_ms = ((double)(t2 - t1)) / 1000.0;

#ifdef COUNT_OSWAPS
    size_t oswap_ap = OSWAP_COUNTER;
#else
    size_t oswap_ap = 0;
#endif

    // ----------------
    // Populate ret
    // ----------------
    if (ret)
    {
        ret->control_bits_time = cb_ms;
        ret->apply_perm_time = ap_ms;
        ret->ptime = cb_ms + ap_ms;

        ret->OSWAP_cb = oswap_cb;
        ret->OSWAP_ap = oswap_ap;
        ret->OSWAP_count = oswap_cb + oswap_ap;
    }

    return cb_ms + ap_ms;
}
