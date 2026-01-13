#include "ButterflyCompaction.hpp"
#include <cstring>
#include <cstdio>
#include <vector>

struct BNPacket
{
    uint32_t value : 1;
    uint32_t tag : 31;
};

static_assert(sizeof(BNPacket) == 4, "BNPacket size must be 4 bytes");

inline void set_bit(uint8_t* arr, size_t idx) {
    arr[idx >> 3] |= (uint8_t)(1u << (idx & 7));
}


static bool PrepareCompactionContext(const bool *selected_list, size_t n, std::vector<uint8_t> &ctrlBits, std::vector<BNPacket> &P, uint32_t &out_Npad, uint32_t &level){
    // 1. 输入检查
    if (!selected_list && n > 0) return false;

    ctrlBits.clear();
    P.clear();

    if (n == 0) return false;

    // 2. Padding
    out_Npad = nextPow2_u32((uint32_t)n);
    if (out_Npad < 2u)
        return false;

    // 3. 计算层数 l
    level = 0;
    for (uint32_t t = out_Npad; t > 1; t >>= 1)
        ++level;

    // 4. 初始化 Packet 数组
    P.resize(out_Npad);
    for (uint32_t i = 0; i < out_Npad; ++i){
        bool is_selected = (i < n) ? selected_list[i] : false;  
        P[i].value = is_selected ? 1u : 0u;
        P[i].tag = 0;
    }

    // 5. 计算 Rank (Tag) - Prefix Sum
    uint32_t rank = 0;
    for (uint32_t i = 0; i < out_Npad; ++i) {
        if (P[i].value) {
            if (rank >= (1u << 31)) return false; // 防止 tag 溢出
            P[i].tag = rank;
            rank++;
        } else {
            P[i].tag = 0;
        }
    }

    // 6. 分配输出控制位空间 (One Byte Per Switch)
    const size_t switchesPerStage = out_Npad / 2;
    const size_t totalSwitches = switchesPerStage * (size_t)level;
    const size_t ctrlBitsLenBytes = (totalSwitches + 7) / 8;

    try {
        // resize 自动初始化为 0
        ctrlBits.resize(ctrlBitsLenBytes, 0); 
    } catch (...) { 
        return false; 
    }

    return true;
}


static void RecGenKernel(std::vector<BNPacket>& P, size_t start, size_t len, uint32_t bit_idx, std::vector<uint8_t> &ctrlBits, size_t& t){

    if (len <= 1) return;

    size_t half = len / 2;

    // 1. 递归左半部分
    RecGenKernel(P, start, half, bit_idx - 1, ctrlBits, t);

    // 2. 递归右半部分
    RecGenKernel(P, start + half, half, bit_idx - 1, ctrlBits, t);

    uint8_t* ctrlBitsPtr = ctrlBits.data();

    // 3. 当前层合并
    for (size_t k = 0; k < half; ++k)
    {
        const size_t i = start + k;
        const size_t j = start + half + k;

        const bool bit_i = ((P[i].tag >> (unsigned)bit_idx) & 1u) != 0;
        const bool bit_j = ((P[j].tag >> (unsigned)bit_idx) & 1u) != 0;

        // c=1 表示交换
        const bool c = (P[i].value && bit_i) || (P[j].value && !bit_j);

        if (c) {
            set_bit(ctrlBitsPtr, t);
        }
        t++;

        oswap_buffer<OSWAP_4>(
            (unsigned char *)&P[i],
            (unsigned char *)&P[j],
            4u,
            (uint8_t)c
        );

    }
}


void RecGenerateControlBits(const bool* selected_list, size_t n, std::vector<uint8_t>& ctrlBits)
{
    FOAV_SAFE_CNTXT(Butterfly_REC_CONTROLBIT, n)

    std::vector<BNPacket> P;
    uint32_t Npad = 0;
    uint32_t level = 0;

    if (!PrepareCompactionContext(selected_list, n, ctrlBits, P, Npad, level)) {
        ctrlBits.clear();
        return;
    }

    size_t t = 0;
    if (level > 0) {
        RecGenKernel(P, 0, Npad, level - 1, ctrlBits, t);
    }
}

static inline void RecApplyButterflyCompact(unsigned char* buf, size_t N, size_t block_size, const std::vector<uint8_t>& ctrlBits)
{
    if (block_size == 4)
    {
        RecApplyCompaction<OSWAP_4>(buf, N, block_size, ctrlBits);
    }
    else if (block_size == 8)
    {
        RecApplyCompaction<OSWAP_8>(buf, N, block_size, ctrlBits);
    }
    else if (block_size == 12)
    {
        RecApplyCompaction<OSWAP_12>(buf, N, block_size, ctrlBits);
    }
    else if (block_size % 16 == 0)
    {
        RecApplyCompaction<OSWAP_16X>(buf, N, block_size, ctrlBits);
    }
    else
    {
        RecApplyCompaction<OSWAP_8_16X>(buf, N, block_size, ctrlBits);
    }
}

static inline void RecApplyTripletCompact(const bool *selected_list, unsigned char* buf, size_t N, size_t block_size)
{
    if (block_size == 4)
    {
        RecTripletCompact<OSWAP_4>(selected_list, buf, N, block_size);
    }
    else if (block_size == 8)
    {
        RecTripletCompact<OSWAP_8>(selected_list, buf, N, block_size);
    }
    else if (block_size == 12)
    {
        RecTripletCompact<OSWAP_12>(selected_list, buf, N, block_size);
    }
    else if (block_size % 16 == 0)
    {
        RecTripletCompact<OSWAP_16X>(selected_list, buf, N, block_size);
    }
    else
    {
        RecTripletCompact<OSWAP_8_16X>(selected_list, buf, N, block_size);
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
    // generateControlBits(selected_list, N, ctrlBits);
    RecGenerateControlBits(selected_list, N, ctrlBits);
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
    // applyButterflyCompact(buffer, N, block_size, ctrlBits);
    RecApplyButterflyCompact(buffer, N, block_size, ctrlBits); 
    //RecApplyTripletCompact(selected_list, buffer, N, block_size);
    ocall_clock(&t2);

    double ap_ms = ((double)(t2 - t1)) / 1000.0;

#ifdef COUNT_OSWAPS
    size_t oswap_ap = OSWAP_COUNTER;
#else
    size_t oswap_ap = 0;
#endif

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