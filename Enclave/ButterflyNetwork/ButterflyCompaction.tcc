#ifndef __BUTTERFLY_COMPACTION_TCC__
#define __BUTTERFLY_COMPACTION_TCC__

#include <cstdint>
#include <cstddef>
#include <vector>
#include <cstring>
#include <cstdio>
#include <cstdlib>

template <OSwap_Style oswap_style>
inline void applyCompaction(unsigned char* buf, size_t n, size_t block_size, const std::vector<uint8_t>& ctrlBits) { 
    FOAV_SAFE2_CNTXT(Butterfly_applyCompaction, n, block_size)

    if (!buf || n < 2) return;

    const uint32_t Npad = nextPow2_u32((uint32_t)n);
    if (Npad < 2u) return;

    uint32_t l = 0;
    {
        uint32_t tmp = Npad;
        while (tmp > 1u) { tmp >>= 1u; ++l; }
    }

    // 1. 计算总开关数
    const size_t totalSwitches = ((size_t)Npad / 2) * (size_t)l;
    
    // 期望大小不再是 Packed Size，而是 totalSwitches
    const size_t expectedSize = totalSwitches;

    if (ctrlBits.size() != expectedSize) {
        printf("applyCompaction: ctrlBits.size=%zu expected=%zu\n",
               ctrlBits.size(), expectedSize);
        return;
    }

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

    size_t t = 0;
    for (uint32_t d = 0; d < l; ++d) {
        const uint32_t stride = (1u << d);

        for (uint32_t b = 0; b < Npad; b += (stride << 1)) {
            for (uint32_t k = 0; k < stride; ++k) {
                const uint32_t i = b + k;
                const uint32_t j = b + k + stride;

                // 【修改点】直接读取 Byte，移除位运算
                const uint8_t c = ctrlBits[t]; 
                t++;

                oswap_buffer<oswap_style>(
                    work + (size_t)block_size * (size_t)i,
                    work + (size_t)block_size * (size_t)j,
                    (uint32_t)block_size,
                    c);
            }
        }
    }

    if (tmp) {
        std::memcpy(buf, tmp, n * block_size);
        free(tmp);
    }
}


template <OSwap_Style oswap_style>
static void RecApplyKernel(unsigned char* work, size_t start, size_t len, size_t block_size, const std::vector<uint8_t>& ctrlBits, size_t& t)
{
    if (len <= 1) return;

    size_t half = len / 2;

    RecApplyKernel<oswap_style>(work, start, half, block_size, ctrlBits, t);
    RecApplyKernel<oswap_style>(work, start + half, half, block_size, ctrlBits, t);

    for (size_t k = 0; k < half; ++k)
    {
        const size_t i = start + k;
        const size_t j = start + half + k;

        const uint8_t c = ctrlBits[t];
        t++;

        oswap_buffer<oswap_style>(
            work + (size_t)block_size * (size_t)i,
            work + (size_t)block_size * (size_t)j,
            (uint32_t)block_size,
            c);
    }
}

template <OSwap_Style oswap_style>
inline void RecApplyCompaction(unsigned char* buf, size_t n, size_t block_size, const std::vector<uint8_t>& ctrlBits) {
    FOAV_SAFE2_CNTXT(Butterfly_RecApply, n, block_size)

    if (!buf || n < 2) return;

    const uint32_t Npad = nextPow2_u32((uint32_t)n);
    if (Npad < 2u) return;

    uint32_t l = 0;
    {
        uint32_t tmp = Npad;
        while (tmp > 1u) { tmp >>= 1u; ++l; }
    }
    const size_t totalSwitches = ((size_t)Npad / 2) * (size_t)l;
    
    // 校验大小为 totalSwitches
    const size_t expectedSize = totalSwitches;

    if (ctrlBits.size() != expectedSize) {
        printf("RecApplyCompaction: ctrlBits.size=%zu expected=%zu\n",
               ctrlBits.size(), expectedSize);
        return;
    }

    unsigned char* work = buf;
    unsigned char* tmp = nullptr;
    if ((uint32_t)n != Npad) {
        tmp = (unsigned char*)malloc((size_t)Npad * block_size);
        if (!tmp) {
            printf("RecApplyCompaction: malloc failed (Npad=%u)\n", Npad);
            return;
        }
        std::memcpy(tmp, buf, n * block_size);
        std::memset(tmp + n * block_size, 0, ((size_t)Npad - n) * block_size);
        work = tmp;
    }

    size_t t = 0;
    RecApplyKernel<oswap_style>(work, 0, Npad, block_size, ctrlBits, t);

    if (tmp) {
        std::memcpy(buf, tmp, n * block_size);
        free(tmp);
    }
}




/**
 * 核心递归内核：RecTripletKernel
 * 对应伪代码: Function RecTriplet(T, start, len, bit)
 * * @param buffer       交错排列的大数组指针
 * @param element_size 单个元素的总大小 (block_size + 8)
 * @param meta_offset  Meta数据起始位置 (等于 block_size)
 * @param bit          当前处理的 bit 位
 */
template <OSwap_Style oswap_style>
static void RecTripletKernel(unsigned char* buffer, size_t start, size_t len, 
                             size_t element_size, size_t meta_offset, uint32_t bit) 
{
    // 1. if len <= 1 then return
    if (len <= 1) return;

    size_t half = len / 2;

    // 3. Step 1: Recurse Depth-First
    RecTripletKernel<oswap_style>(buffer, start, half, element_size, meta_offset, bit - 1);
    RecTripletKernel<oswap_style>(buffer, start + half, half, element_size, meta_offset, bit - 1);

    // 6. Step 2: Merge (Single Swap per Butterfly)
    for (size_t k = 0; k < half; ++k) {
        const size_t i = start + k;
        const size_t j = start + half + k;

        // 计算物理地址
        unsigned char* ptr_i = buffer + i * element_size;
        unsigned char* ptr_j = buffer + j * element_size;

        // 获取 Meta 数据 (Tag 和 Act)
        // 内存布局: [ Data ... | Tag (4B) | Act (4B) ]
        // meta_offset 指向 Tag 的位置
        uint32_t* meta_i_ptr = (uint32_t*)(ptr_i + meta_offset);
        uint32_t* meta_j_ptr = (uint32_t*)(ptr_j + meta_offset);

        // 读取 Tag (偏移 0) 和 Act (偏移 1, 即 +4 bytes)
        const uint32_t tag_i = meta_i_ptr[0]; 
        const uint32_t act_i = meta_i_ptr[1]; 
        const uint32_t tag_j = meta_j_ptr[0]; 
        const uint32_t act_j = meta_j_ptr[1]; 

        // 10. Compute control logic
        const uint32_t bit_i = (tag_i >> bit) & 1u;
        const uint32_t bit_j = (tag_j >> bit) & 1u;

        const uint32_t down_i = act_i & bit_i;
        const uint32_t up_j   = act_j & (bit_j ^ 1u);
        const uint8_t  c      = (uint8_t)((down_i | up_j) & 1u);

        // 16. Execute Single OSwap (Moves Data, Tag, Act together)
        // 交换整个元素块 (Data + Meta)
        oswap_buffer<oswap_style>(
            ptr_i,
            ptr_j,
            (uint32_t)element_size,
            c
        );
    }
}

/**
 * 三元组压缩入口函数：RecTripletCompact
 * 对应伪代码: REC_TRIPLET_COMPACT(A, D)
 */
template <OSwap_Style oswap_style>
inline void RecTripletCompact(const bool *selected_list, unsigned char* buf, size_t n, size_t block_size) {
    FOAV_SAFE2_CNTXT(Butterfly_RecTriplet, n, block_size)

    if (!buf || n < 1) return;

    // 2. l <- ceil(log2 n), N <- 2^l
    const uint32_t N = nextPow2_u32((uint32_t)n);
    uint32_t l = 0;
    {
        uint32_t tmp = N;
        while (tmp > 1u) { tmp >>= 1u; ++l; }
    }

    // -----------------------------------------------------------
    // Memory Layout Strategy:
    // [ Payload (block_size) ] [ Tag (4 bytes) ] [ Act (4 bytes) ]
    // Total overhead per element = 8 bytes
    // -----------------------------------------------------------
    
    const size_t tag_size = 4;
    const size_t act_size = 4;
    const size_t meta_size = tag_size + act_size; // 8 bytes
    const size_t element_size = block_size + meta_size;
    
    // 3. Create array T of length N
    // 分配大缓冲区存放交错数据
    unsigned char* T = (unsigned char*)malloc(N * element_size);
    if (!T) {
        printf("RecTripletCompact: malloc failed (N=%u)\n", N);
        return;
    }

    // 4. Phase 0: Initialization & Packing
    for (size_t i = 0; i < N; ++i) {
        unsigned char* elem_ptr = T + i * element_size;
        
        // 指向 Meta 区域
        uint32_t* meta_ptr = (uint32_t*)(elem_ptr + block_size);
        // meta_ptr[0] is Tag, meta_ptr[1] is Act

        if (i < n) {
            // 6. Copy Data: D[i]
            std::memcpy(elem_ptr, buf + i * block_size, block_size);
            // Set Meta: tag=0, act=A[i]
            meta_ptr[0] = 0;
            meta_ptr[1] = selected_list[i] ? 1 : 0;
        } else {
            // 8. Padding
            std::memset(elem_ptr, 0, block_size);
            meta_ptr[0] = 0;
            meta_ptr[1] = 0;
        }
    }

    // 10. Phase 1: Oblivious Ranking
    uint32_t rank = 0;
    for (size_t i = 0; i < N; ++i) {
        unsigned char* elem_ptr = T + i * element_size;
        uint32_t* meta_ptr = (uint32_t*)(elem_ptr + block_size);
        
        // 引用以便读写
        uint32_t& tag = meta_ptr[0];
        uint32_t& act = meta_ptr[1];

        // 11. T[i].tag <- rank * T[i].act
        tag = rank * act;

        // 12. rank <- rank + T[i].act
        rank += act;
    }

    // 13. Call RecTriplet (Phase 2: Unified Permutation)
    if (l > 0) {
        // bit 初始为 l - 1
        // 传入 block_size 作为 meta 在元素内的偏移量
        RecTripletKernel<oswap_style>(T, 0, N, element_size, block_size, l - 1);
    }

    // 14. Unpack result
    // 将 T 中的数据部分拷回 buf
    // 只需要拷贝前 n 个块 (实际上只有前 rank 个是有效的，但为了保持接口行为一致，拷贝 n 个)
    for (size_t i = 0; i < n; ++i) {
        unsigned char* src_ptr = T + i * element_size; // 只需要由 T 指向数据段开头
        std::memcpy(buf + i * block_size, src_ptr, block_size);
    }

    free(T);
}


#endif // __BUTTERFLY_COMPACTION_TCC__