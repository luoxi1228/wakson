#ifndef __BUTTERFLY_COMPACTION_TCC__
#define __BUTTERFLY_COMPACTION_TCC__

#include <cstdint>
#include <cstddef>
#include <vector>
#include <cstring>
#include <cstdio>
#include <cstdlib>


// 读取压缩位
inline bool get_bit(const std::vector<uint8_t>& bits, size_t idx) {
    return (bits[idx >> 3] >> (idx & 7)) & 1u;
}

// 专门用于快速读取位的流对象
struct BitStream {
    const uint8_t* data;  // 指向 vector 数据的裸指针
    size_t byte_idx;      // 当前读取到的字节索引
    int bit_idx;          // 当前字节内的位偏移 (0-7)
    uint8_t current_byte; // 缓存当前字节，减少内存访问

    // 构造函数
    BitStream(const std::vector<uint8_t>& bits) {
        data = bits.data();
        byte_idx = 0;
        bit_idx = 0;
        // 预加载第一个字节（假设 bits 不为空）
        if (!bits.empty()) {
            current_byte = data[0];
        }
    }

    // 极速读取下一个位 (Force Inline)
    inline bool next() {
        // 直接从寄存器缓存中提取位
        bool val = (current_byte >> bit_idx) & 1u;
        
        bit_idx++;
        // 分支预测非常友好，每8次才触发一次
        if (bit_idx == 8) {
            byte_idx++;
            current_byte = data[byte_idx]; // 只有这里才会访问内存
            bit_idx = 0;
        }
        return val;
    }
};

template <OSwap_Style oswap_style>
static void RecApplyKernel(
    unsigned char* work,
    size_t start,
    size_t len,
    size_t block_size,
    const std::vector<uint8_t>& ctrlBits,
    BitStream& bs)
{
    if (len <= 1) return;

    const size_t half = len / 2;

    RecApplyKernel<oswap_style>(work, start, half, block_size, ctrlBits, bs);
    RecApplyKernel<oswap_style>(work, start + half, half, block_size, ctrlBits, bs);

    for (size_t k = 0; k < half; ++k)
    {
        const size_t i = start + k;
        const size_t j = start + half + k;

        bool c = bs.next();

        oswap_buffer<oswap_style>(
            work + block_size * i,
            work + block_size * j,
            (uint32_t)block_size,
            (uint8_t)c);

    }
}


template <OSwap_Style oswap_style>
inline void RecApplyCompaction(
    unsigned char* buf,
    size_t n,
    size_t block_size,
    const std::vector<uint8_t>& ctrlBits)
{
    FOAV_SAFE2_CNTXT(Butterfly_RecApply, n, block_size)
    if (!buf || n < 2) return;

    const uint32_t Npad = nextPow2_u32((uint32_t)n);
    if (Npad < 2u) return;

    uint32_t l = 0;
    for (uint32_t tmp = Npad; tmp > 1u; tmp >>= 1u) ++l;

    const size_t totalSwitches = (size_t)(Npad / 2) * (size_t)l;
    const size_t expectedSize = (totalSwitches + 7) / 8;
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

    BitStream bs(ctrlBits);

    RecApplyKernel<oswap_style>(work, 0, Npad, block_size, ctrlBits, bs);


    if (tmp) {
        std::memcpy(buf, tmp, n * block_size);
        free(tmp);
    }
}



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