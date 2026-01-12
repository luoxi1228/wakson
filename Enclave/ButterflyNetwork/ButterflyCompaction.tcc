#ifndef __BUTTERFLY_COMPACTION_TCC__
#define __BUTTERFLY_COMPACTION_TCC__

#include <cstdint>
#include <cstddef>
#include <vector>
#include <cstring>
#include <cstdio>
#include <cstdlib>

template <OSwap_Style oswap_style>
inline void applyCompaction(unsigned char* buf,
                            size_t n,
                            size_t block_size,
                            const std::vector<uint8_t>& ctrlBits) { // 注意：这里的 vector 存储的是压缩后的数据
    FOAV_SAFE2_CNTXT(Butterfly_applyCompaction, n, block_size)

    if (!buf || n < 2) return;

    const uint32_t Npad = nextPow2_u32((uint32_t)n);
    if (Npad < 2u) return;

    // l = log2(Npad)
    uint32_t l = 0;
    {
        uint32_t tmp = Npad;
        while (tmp > 1u) { tmp >>= 1u; ++l; }
    }

    // 1. 计算总开关数 (Total Switches)
    const size_t totalSwitches = ((size_t)Npad / 2) * (size_t)l;
    
    // 2. 【修改点一】计算期望的字节数 (Packed Size)
    // (totalSwitches + 7) / 8 等价于位运算 (total + 7) >> 3
    const size_t expectedPackedSize = (totalSwitches + 7) >> 3;

    if (ctrlBits.size() != expectedPackedSize) {
        printf("applyCompaction: ctrlBits.size=%zu expectedPacked=%zu (totalSwitches=%zu)\n",
               ctrlBits.size(), expectedPackedSize, totalSwitches);
        return;
    }

    // 如果 n 不是 2^k，必须 pad，否则 online 阶段会越界访问 i/j
    unsigned char* work = buf;
    unsigned char* tmp = nullptr;
    if ((uint32_t)n != Npad) {
        tmp = (unsigned char*)malloc((size_t)Npad * block_size);
        if (!tmp) {
            printf("applyCompaction: malloc failed (Npad=%u)\n", Npad);
            return;
        }
        std::memcpy(tmp, buf, n * block_size);
        std::memset(tmp + n * block_size, 0, ((size_t)Npad - n) * block_size);
        work = tmp;
    }

    // online: 按 ctrlBits 的线性顺序应用
    size_t t = 0;
    for (uint32_t d = 0; d < l; ++d) {
        const uint32_t stride = (1u << d);

        for (uint32_t b = 0; b < Npad; b += (stride << 1)) {
            for (uint32_t k = 0; k < stride; ++k) {
                const uint32_t i = b + k;
                const uint32_t j = b + k + stride;

                // 3. 【修改点二】位解压读取 (Bit Unpacking)
                // 从压缩数组中提取第 t 个 bit
                // t >> 3 确定字节索引，t & 7 确定位偏移
                const uint8_t c = (uint8_t)((ctrlBits[t >> 3] >> (t & 7)) & 1u);
                t++;

                oswap_buffer<oswap_style>(
                    work + (size_t)block_size * (size_t)i,
                    work + (size_t)block_size * (size_t)j,
                    (uint32_t)block_size,
                    c);
            }
        }
    }

    // copy back if padded
    if (tmp) {
        std::memcpy(buf, tmp, n * block_size);
        free(tmp);
    }
}


/**
 * 递归核心内核
 * @param t: 全局 bit 计数器引用，确保按顺序消耗压缩的控制位
 */
template <OSwap_Style oswap_style>
static void RecApplyKernel(unsigned char* work,
                           size_t start,
                           size_t len,
                           size_t block_size,
                           const std::vector<uint8_t>& ctrlBits,
                           size_t& t)
{
    // Base Case: 长度为 1 或 0 时停止
    if (len <= 1) return;

    size_t half = len / 2;

    // 1. 递归处理左半部分 (Depth-First)
    RecApplyKernel<oswap_style>(work, start, half, block_size, ctrlBits, t);

    // 2. 递归处理右半部分
    RecApplyKernel<oswap_style>(work, start + half, half, block_size, ctrlBits, t);

    // 3. 当前层合并 (Merge)
    // 对应递归生成中的 Merge 阶段 (处理 stride = half 的层)
    for (size_t k = 0; k < half; ++k)
    {
        const size_t i = start + k;
        const size_t j = start + half + k;

        // 【位解压读取】
        // 逻辑与迭代版完全一致：(arr[byte_idx] >> bit_offset) & 1
        const uint8_t c = (uint8_t)((ctrlBits[t >> 3] >> (t & 7)) & 1u);
        t++;

        // 执行不经意交换
        oswap_buffer<oswap_style>(
            work + (size_t)block_size * (size_t)i,
            work + (size_t)block_size * (size_t)j,
            (uint32_t)block_size,
            c);
    }
}

template <OSwap_Style oswap_style>
inline void RecApplyCompaction(unsigned char* buf,
                               size_t n,
                               size_t block_size,
                               const std::vector<uint8_t>& ctrlBits) 
{
    FOAV_SAFE2_CNTXT(Butterfly_RecApply, n, block_size)

    if (!buf || n < 2) return;

    const uint32_t Npad = nextPow2_u32((uint32_t)n);
    if (Npad < 2u) return;

    // 1. 校验 Control Bits 大小 (与迭代版逻辑一致)
    // l = log2(Npad)
    uint32_t l = 0;
    {
        uint32_t tmp = Npad;
        while (tmp > 1u) { tmp >>= 1u; ++l; }
    }
    const size_t totalSwitches = ((size_t)Npad / 2) * (size_t)l;
    const size_t expectedPackedSize = (totalSwitches + 7) >> 3;

    if (ctrlBits.size() != expectedPackedSize) {
        printf("RecApplyCompaction: ctrlBits.size=%zu expectedPacked=%zu\n",
               ctrlBits.size(), expectedPackedSize);
        return;
    }

    // 2. Padding 处理
    // 递归算法依赖 strict power-of-2 分割，必须使用 Padding 避免越界
    unsigned char* work = buf;
    unsigned char* tmp = nullptr;
    if ((uint32_t)n != Npad) {
        tmp = (unsigned char*)malloc((size_t)Npad * block_size);
        if (!tmp) {
            printf("RecApplyCompaction: malloc failed (Npad=%u)\n", Npad);
            return;
        }
        // 拷贝原数据
        std::memcpy(tmp, buf, n * block_size);
        // Padding 部分清零
        std::memset(tmp + n * block_size, 0, ((size_t)Npad - n) * block_size);
        work = tmp;
    }

    // 3. 执行递归内核
    // 这里的 t 是全局计数器，会在递归过程中被引用修改
    size_t t = 0;
    
    // 从整个数组范围 (0 到 Npad) 开始递归
    RecApplyKernel<oswap_style>(work, 0, Npad, block_size, ctrlBits, t);

    // 4. 写回数据 (如果用了 Padding)
    if (tmp) {
        std::memcpy(buf, tmp, n * block_size);
        free(tmp);
    }
}




#endif // __BUTTERFLY_COMPACTION_TCC__