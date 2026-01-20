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

// 迭代版
void generateControlBits(const bool *selected_list, size_t n, std::vector<uint8_t> &ctrlBits)
{
    FOAV_SAFE_CNTXT(Butterfly_CONTROLBIT, n)

    std::vector<BNPacket> P;
    uint32_t Npad = 0;
    uint32_t l = 0;

    if (!PrepareCompactionContext(selected_list, n, ctrlBits, P, Npad, l)) {
        return;
    }

    size_t t = 0;                 // 全局开关计数器
    uint8_t* ctrlBitsPtr = ctrlBits.data();

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

                const uint8_t c = (uint8_t)((P[i].value & bit_i) |
                                            (P[j].value & (bit_j ^ 1u)));

                if (c) set_bit(ctrlBitsPtr, t);
                ++t;

                oswap_buffer<OSWAP_4>(
                    (unsigned char *)&P[i],
                    (unsigned char *)&P[j],
                    4u,
                    c);
            }
        }
    }
}
static inline void applyButterflyCompact(unsigned char *buf, size_t N, size_t block_size, const std::vector<uint8_t> &ctrlBits)
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


// 递归版
static void generateControlBitsRec_inner(std::vector<BNPacket>& P, size_t start, size_t len, uint32_t bit_idx, std::vector<uint8_t> &ctrlBits, size_t& t){

    if (len <= 1) return;

    size_t half = len / 2;

    // 1. 递归左半部分
    generateControlBitsRec_inner(P, start, half, bit_idx - 1, ctrlBits, t);

    // 2. 递归右半部分
    generateControlBitsRec_inner(P, start + half, half, bit_idx - 1, ctrlBits, t);

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
void generateControlBitsRec(const bool* selected_list, size_t n, std::vector<uint8_t>& ctrlBits)
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
        generateControlBitsRec_inner(P, 0, Npad, level - 1, ctrlBits, t);
    }
}
static inline void applyButterflyCompactRec(unsigned char* buf, size_t N, size_t block_size, const std::vector<uint8_t>& ctrlBits)
{
    if (block_size == 4)
    {
        applyCompactionRec<OSWAP_4>(buf, N, block_size, ctrlBits);
    }
    else if (block_size == 8)
    {
        applyCompactionRec<OSWAP_8>(buf, N, block_size, ctrlBits);
    }
    else if (block_size == 12)
    {
        applyCompactionRec<OSWAP_12>(buf, N, block_size, ctrlBits);
    }
    else if (block_size % 16 == 0)
    {
        applyCompactionRec<OSWAP_16X>(buf, N, block_size, ctrlBits);
    }
    else
    {
        applyCompactionRec<OSWAP_8_16X>(buf, N, block_size, ctrlBits);
    }
}
static inline void applyButterflyCompactRecIter( unsigned char* buf, size_t N, size_t block_size, const std::vector<uint8_t>& ctrlBits)
{
    if (block_size == 4) {
        applyCompactionRecIter<OSWAP_4>(buf, N, block_size, ctrlBits);
    } else if (block_size == 8) {
        applyCompactionRecIter<OSWAP_8>(buf, N, block_size, ctrlBits);
    } else if (block_size == 12) {
        applyCompactionRecIter<OSWAP_12>(buf, N, block_size, ctrlBits);
    } else if (block_size % 16 == 0) {
        applyCompactionRecIter<OSWAP_16X>(buf, N, block_size, ctrlBits);
    } else {
        applyCompactionRecIter<OSWAP_8_16X>(buf, N, block_size, ctrlBits);
    }
}


//OR 版
static void generateControlBitsOR_inner2power(const bool* M, size_t start, size_t n, uint32_t z, std::vector<uint8_t>& Cbits, size_t& t){
    if (n <= 1) return;

    // m = Sum(M[0..n/2-1])  (左半部分 sum)
    const size_t half = n >> 1;
    const uint32_t m = range_sum_bool(M, start, start + half);

    if (n == 2) {
        // bit = ((1 - M0) & M1) ^ z
        const bool M0 = M[start];
        const bool M1 = M[start + 1];
        const bool bit = (((!M0) && M1) ? 1 : 0) ^ (z & 1u);
        append_bit(Cbits, t, bit);
        return;
    }

    // 递归（顺序必须与 Online 一致）
    generateControlBitsOR_inner2power(M, start,         half, (uint32_t)(z % half),              Cbits, t);
    generateControlBitsOR_inner2power(M, start + half,  half, (uint32_t)((z + m) % half),        Cbits, t);

    // s = (((z mod half) + m) >= half) ^ (z >= half)
    const bool s = ((((uint32_t)(z % half) + m) >= (uint32_t)half) ? 1 : 0) ^ ((z >= (uint32_t)half) ? 1 : 0);

    const uint32_t threshold = (uint32_t)((z + m) % half);
    for (uint32_t i = 0; i < (uint32_t)half; ++i) {
        const bool bit = s ^ (i >= threshold);
        append_bit(Cbits, t, bit);
    }
}
static void generateControlBitsOR_inner(const bool* M, size_t start, size_t n, std::vector<uint8_t>& Cbits, size_t& t){
    if (n == 0) return;

    // n1 = largest power of two <= n
    // n2 = n - n1
    size_t n1 = 1;
    while ((n1 << 1) <= n) n1 <<= 1;
    const size_t n2 = n - n1;

    // m = Sum(M[0..n2-1]) (prefix sum over the non-power-of-two part)
    const uint32_t m = range_sum_bool(M, start, start + n2);

    // 递归处理 prefix (n2)
    generateControlBitsOR_inner(M, start, n2, Cbits, t);

    // 处理 suffix (n1)，offset z = n1 - n2 + m
    const uint32_t z = (uint32_t)(n1 - n2) + m;
    generateControlBitsOR_inner2power(M, start + n2, n1, z, Cbits, t);

    // 最终合并 loop：for i in [0..n2-1] append (i > m)
    for (uint32_t i = 0; i < (uint32_t)n2; ++i) {
        const bool bit = (i > m);
        append_bit(Cbits, t, bit);
    }
}
void generateControlBitsOR(const bool* selected_list, size_t n, std::vector<uint8_t>& ctrlBits)
{
    FOAV_SAFE_CNTXT(OR_COMPACT_CONTROLBIT, n);
    
    ctrlBits.clear();
    
    if (!selected_list || n == 0) return;

    const size_t roughBits = (n * 2); // 粗略：先给少点也没关系，append_bit 会自动扩容

    ctrlBits.reserve((roughBits + 7) / 8);

    size_t t = 0;
    generateControlBitsOR_inner(selected_list, 0, n, ctrlBits, t);

    // 末尾把 vector 修到刚好需要的 byte（append_bit 已经保证足够）
    const size_t needBytes = (t + 7) / 8;
    ctrlBits.resize(needBytes, 0);

}
static inline void applyButterflyCompactOR( unsigned char* buf, size_t N, size_t block_size, const std::vector<uint8_t>& ctrlBits)
{
    if (block_size == 4) {
        applyCompactionOR<OSWAP_4>(buf, N, block_size, ctrlBits);
    } else if (block_size == 8) {
        applyCompactionOR<OSWAP_8>(buf, N, block_size, ctrlBits);
    } else if (block_size == 12) {
        applyCompactionOR<OSWAP_12>(buf, N, block_size, ctrlBits);
    } else if (block_size % 16 == 0) {
        applyCompactionOR<OSWAP_16X>(buf, N, block_size, ctrlBits);
    } else {
        applyCompactionOR<OSWAP_8_16X>(buf, N, block_size, ctrlBits);
    }
}


//递归合并版
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



static void OptGenerateControlBitsOR_inner2power( const bool* M, const uint32_t* ps, uint32_t start, uint32_t n, uint32_t z, uint8_t* out, size_t& t)
{
    if (n <= 1) return;

    uint32_t half = n >> 1;
    uint32_t m = ps[start + half] - ps[start]; // sum(left half)

    if (n == 2) {
        bool M0 = M[start];
        bool M1 = M[start + 1];
        uint8_t bit = (uint8_t)(((!M0 && M1) ? 1 : 0) ^ (z & 1u));
        out[t++] = bit;
        return;
    }

    // recurse left/right in the same order as apply
    OptGenerateControlBitsOR_inner2power(M, ps, start,        half, (z & (half - 1u)), out, t);
    OptGenerateControlBitsOR_inner2power(M, ps, start + half, half, ((z + m) & (half - 1u)), out, t);

    // merge bits: bit(i) = s ^ (i >= threshold)
    uint32_t zmod = (z & (half - 1u));
    uint8_t offset_right = (uint8_t)(z >= half);
    uint8_t left_wrapped = (uint8_t)((zmod + m) >= half);
    uint8_t s = (uint8_t)(left_wrapped ^ offset_right);
    uint32_t threshold = (uint32_t)((z + m) & (half - 1u));

    // fill [0, threshold) with s; [threshold, half) with s^1
    if (threshold) {
        std::memset(out + t, s, threshold);
        t += threshold;
    }
    uint32_t rest = half - threshold;
    if (rest) {
        std::memset(out + t, (uint8_t)(s ^ 1u), rest);
        t += rest;
    }
}
static void OptGenerateControlBitsOR_inner( const bool* M, const uint32_t* ps, uint32_t start, uint32_t n, uint8_t* out, size_t& t)
{
    if (n <= 1) return;

    if (n == 2) {
        // 对应 offset=0 的 2power base case
        bool M0 = M[start];
        bool M1 = M[start + 1];
        out[t++] = (uint8_t)((!M0 && M1) ? 1 : 0);
        return;
    }

    if (is_pow2_u32(n)) {
        OptGenerateControlBitsOR_inner2power(M, ps, start, n, 0u, out, t);
        return;
    }

    // general split: n1=pow2<=n, n2=n-n1
    uint32_t n1 = pow2_le_u32(n);
    uint32_t n2 = n - n1;

    // m = sum(prefix n2)
    uint32_t m = ps[start + n2] - ps[start];

    // 1) prefix general
    OptGenerateControlBitsOR_inner(M, ps, start, n2, out, t);

    // 2) suffix 2power with offset z = (n1 - n2) + m
    uint32_t z = (n1 - n2) + m;
    OptGenerateControlBitsOR_inner2power(M, ps, start + n2, n1, z, out, t);

    // 3) final merge bits for i in [0..n2-1]
    // original: bit = (i > m)
    // => zeros length = min(n2, m+1), ones afterwards
    uint32_t zeros = n2;
    uint32_t m1 = m + 1u; // careful overflow not possible here
    if (m1 < zeros) zeros = m1;

    if (zeros) {
        std::memset(out + t, 0u, zeros);
        t += zeros;
    }
    uint32_t ones = n2 - zeros;
    if (ones) {
        std::memset(out + t, 1u, ones);
        t += ones;
    }
}
void OptGenerateControlBitsOR(const bool* selected_list, size_t n, std::vector<uint8_t>& ctrlBits)
{
    FOAV_SAFE_CNTXT(OR_COMPACT_CONTROLBIT, n);

    ctrlBits.clear();
    if (!selected_list || n == 0) return;

    uint32_t N = (uint32_t)n;
    if (N < 2u) return;

    //计算总开关数
    size_t totalSwitches = count_or_switches(N);
    ctrlBits.resize(totalSwitches);

    // prefix sum for O(1) range sum
    std::vector<uint32_t> ps(N + 1);
    ps[0] = 0;
    for (uint32_t i = 0; i < N; ++i) {
        ps[i + 1] = ps[i] + (uint32_t)selected_list[i];
    }

    size_t t = 0;
    OptGenerateControlBitsOR_inner(selected_list, ps.data(), 0u, N, ctrlBits.data(), t);
}
static inline void OptApplyButterflyCompactOR( unsigned char* buf, size_t N, size_t block_size, const std::vector<uint8_t>& ctrlBits)
{
    if (block_size == 4) {
        OptApplyCompactionOR<OSWAP_4>(buf, N, block_size, ctrlBits);
    } else if (block_size == 8) {
        OptApplyCompactionOR<OSWAP_8>(buf, N, block_size, ctrlBits);
    } else if (block_size == 12) {
        OptApplyCompactionOR<OSWAP_12>(buf, N, block_size, ctrlBits);
    } else if (block_size % 16 == 0) {
        OptApplyCompactionOR<OSWAP_16X>(buf, N, block_size, ctrlBits);
    } else {
        OptApplyCompactionOR<OSWAP_8_16X>(buf, N, block_size, ctrlBits);
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
    //generateControlBits(selected_list, N, ctrlBits);
    //generateControlBitsRec(selected_list, N, ctrlBits);
    //generateControlBitsOR(selected_list, N, ctrlBits);
    OptGenerateControlBitsOR(selected_list, N, ctrlBits);
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
    //applyButterflyCompactRecIter(buffer, N, block_size, ctrlBits);
    //applyButterflyCompactRec(buffer, N, block_size, ctrlBits); 
    //applyButterflyCompactOR(buffer, N, block_size, ctrlBits);
    //RecApplyTripletCompact(selected_list, buffer, N, block_size);
    OptApplyButterflyCompactOR(buffer, N, block_size, ctrlBits);
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