#include "ButterflyCompaction.hpp"
#include <cstring>
#include <cstdio>
#include <vector>

struct BNPacket
{
    uint32_t value;
    uint32_t tag;
};

static_assert(sizeof(BNPacket) == 8, "BNPacket must be 8 bytes");

// 辅助宏：向位压缩数组写入第 i 位 (设为 1)
#define SET_BIT(arr, i) ((arr)[(i) >> 3] |= (uint8_t)(1u << ((i) & 7)))


/**
 * 负责通用的初始化工作：
 * - 初始化 P 数组
 * - 计算 Rank
 * - 分配位压缩的 ctrlBits 空间
 */
static bool PrepareCompactionContext(const bool *selected_list, size_t n, std::vector<uint8_t> &ctrlBits, std::vector<BNPacket> &P, uint32_t &out_Npad, uint32_t &out_l){
    // 1. 输入检查
    ctrlBits.clear();
    if (!selected_list || n == 0)
        return false;

    // 2. 计算 Padding
    out_Npad = nextPow2_u32((uint32_t)n);
    if (out_Npad < 2u)
        return false;

    // 3. 计算层数 l
    out_l = 0;
    for (uint32_t t = out_Npad; t > 1; t >>= 1)
        ++out_l;

    // 4. 初始化 Packet 数组
    P.resize(out_Npad);
    for (uint32_t i = 0; i < out_Npad; ++i)
    {
        // 【保持原逻辑】直接读取 bool 数组
        // 如果 i < n，读取 selected_list[i]；Padding 部分补 0
        uint32_t is_selected = (i < n && selected_list[i]) ? 1u : 0u;
        
        P[i].value = is_selected;
        P[i].tag = 0;
    }

    // 5. 计算 Rank (Tag) - Prefix Sum
    uint32_t rank = 0;
    for (uint32_t i = 0; i < out_Npad; ++i)
    {
        const uint32_t v = P[i].value;
        P[i].tag = rank * v; // value为0时tag为0，value为1时tag为rank
        rank += v;
    }

    // 6. 分配输出控制位空间 (Bit-Packed)
    // 总开关数 = (N/2) * l
    // 字节数 = (TotalBits + 7) / 8
    const size_t switchesPerStage = out_Npad / 2;
    const size_t totalSwitches = switchesPerStage * out_l;
    const size_t packedSize = (totalSwitches + 7) >> 3; 
    
    ctrlBits.assign(packedSize, 0);

    return true;
}

static void GenKernel( std::vector<BNPacket>& P, uint32_t Npad, uint32_t l, std::vector<uint8_t> &ctrlBits){
    size_t t = 0; // 全局开关计数器
    for (uint32_t d = 0; d < l; ++d)
    {
        const uint32_t stride = (1u << d);
        // 遍历所有 Butterfly Block
        for (uint32_t b = 0; b < Npad; b += (stride << 1))
        {
            // Block 内部遍历
            for (uint32_t k = 0; k < stride; ++k)
            {
                const uint32_t i = b + k;
                const uint32_t j = b + k + stride;

                const uint32_t bit_i = (P[i].tag >> d) & 1u;
                const uint32_t bit_j = (P[j].tag >> d) & 1u;

                const uint8_t c = (uint8_t)((P[i].value & bit_i) |
                                            (P[j].value & (bit_j ^ 1u)));

                // 【位压缩写入】
                if (c) {
                    SET_BIT(ctrlBits, t);
                }
                t++;

                oswap_buffer<OSWAP_8>(
                    (unsigned char *)&P[i],
                    (unsigned char *)&P[j],
                    8u,
                    c);
            }
        }
    }
}

static void RecGenKernel(std::vector<BNPacket>& P, size_t start, size_t len, int bit_idx, std::vector<uint8_t>& ctrlBits, size_t& t){
    // Base Case
    if (len <= 1) return;
    
    size_t half = len / 2;

    // 1. 递归左半部分 (bit_idx - 1)
    RecGenKernel(P, start, half, bit_idx - 1, ctrlBits, t);

    // 2. 递归右半部分 (bit_idx - 1)
    RecGenKernel(P, start + half, half, bit_idx - 1, ctrlBits, t);

    // 3. 当前层合并 (Merge)
    // 此时 bit 0 ... bit_idx-1 已经处理完毕，现在根据 bit_idx 进行交换
    for (size_t k = 0; k < half; ++k)
    {
        const size_t i = start + k;
        const size_t j = start + half + k;

        uint32_t d = (uint32_t)bit_idx;

        const uint32_t bit_i = (P[i].tag >> d) & 1u;
        const uint32_t bit_j = (P[j].tag >> d) & 1u;

        const uint8_t c = (uint8_t)((P[i].value & bit_i) |
                                    (P[j].value & (bit_j ^ 1u)));

        // 【位压缩写入】
        if (c) {
            SET_BIT(ctrlBits, t);
        }
        t++;

        oswap_buffer<OSWAP_8>(
            (unsigned char *)&P[i],
            (unsigned char *)&P[j],
            8u,
            c);
    }
}


// 迭代版接口
void generateControlBits(const bool *selected_list, size_t n, std::vector<uint8_t> &ctrlBits)
{
    FOAV_SAFE_CNTXT(Butterfly_CONTROLBIT, n)

    std::vector<BNPacket> P;
    uint32_t Npad = 0;
    uint32_t l = 0;

    // 调用公共准备函数
    if (!PrepareCompactionContext(selected_list, n, ctrlBits, P, Npad, l)) {
        return;
    }

    // 执行迭代逻辑
    GenKernel(P, Npad, l, ctrlBits);
}
// 递归版接口
void RecGenerateControlBits(const bool* selected_list, size_t n, std::vector<uint8_t>& ctrlBits)
{
    FOAV_SAFE_CNTXT(Butterfly_REC_CONTROLBIT, n)

    std::vector<BNPacket> P;
    uint32_t Npad = 0;
    uint32_t l = 0;

    // 调用公共准备函数
    if (!PrepareCompactionContext(selected_list, n, ctrlBits, P, Npad, l)) {
        return;
    }

    // 执行递归逻辑
    size_t t = 0;
    if (l > 0) {
        // 从最高层 l-1 开始递归，对应整个数组
        RecGenKernel(P, 0, Npad, (int)l - 1, ctrlBits, t);
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



