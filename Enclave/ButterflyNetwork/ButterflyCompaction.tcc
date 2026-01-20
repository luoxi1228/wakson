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
    const uint8_t* data;
    size_t nbytes;
    size_t byte_idx;
    int bit_idx;
    uint8_t current_byte;

    BitStream(const std::vector<uint8_t>& bits)
        : data(bits.data()),
          nbytes(bits.size()),
          byte_idx(0),
          bit_idx(0),
          current_byte((nbytes > 0) ? bits[0] : 0u) {}

    inline uint8_t next() {
        uint8_t val = (uint8_t)((current_byte >> bit_idx) & 1u);

        ++bit_idx;
        if (bit_idx == 8) {
            bit_idx = 0;
            ++byte_idx;
            current_byte = (byte_idx < nbytes) ? data[byte_idx] : 0u; 
        }
        return val;
    }
};

static inline size_t pow2_lt_size(size_t n) {
    size_t p = 1;
    while ((p << 1) <= n) p <<= 1;
    return p;
}

//迭代
template <OSwap_Style oswap_style>
inline void applyCompaction(unsigned char* buf, size_t n, size_t block_size, const std::vector<uint8_t>& ctrlBits)
{
    FOAV_SAFE2_CNTXT(Butterfly_applyCompaction, n, block_size)
    if (!buf || n < 2) return;

    const uint32_t Npad = nextPow2_u32((uint32_t)n);
    if (Npad < 2u) return;

    uint32_t l = 0;
    for (uint32_t tmpN = Npad; tmpN > 1u; tmpN >>= 1u) ++l;

    const size_t totalSwitches = ((size_t)Npad / 2) * (size_t)l;
    const size_t expectedSize  = (totalSwitches + 7) / 8; 
    if (ctrlBits.size() != expectedSize) {
        printf("applyCompaction: ctrlBits.size=%zu expected=%zu\n",
               ctrlBits.size(), expectedSize);
        return;
    }

    unsigned char* work = buf;
    unsigned char* tmp  = nullptr;
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

    BitStream bs(ctrlBits);

    //const uint8_t* ctrl = ctrlBits.data(); // 裸指针读控制位更快
    //size_t t = 0;


    for (uint32_t d = 0; d < l; ++d) {
        const uint32_t stride = (1u << d);
        const size_t stride_bytes = (size_t)stride * block_size;
        const uint32_t step = stride << 1;

        for (uint32_t b = 0; b < Npad; b += step) {
            unsigned char* base = work + (size_t)b * block_size;

            unsigned char* ptr_i = base;
            unsigned char* ptr_j = base + stride_bytes;

            // k: 0..stride-1
            for (uint32_t k = 0; k < stride; ++k) {
                bool c = bs.next();
                //const uint8_t c = (uint8_t)(ctrl[t++] & 1u); // 保证是 0/1
                oswap_buffer<oswap_style>(ptr_i, ptr_j, block_size, c);
                ptr_i += block_size;
                ptr_j += block_size;
            }
        }
    }

    if (tmp) {
        std::memcpy(buf, tmp, n * block_size);
        free(tmp);
    }
}


// 递归
template <OSwap_Style oswap_style>
static void applyCompactionRec_inner(unsigned char* current_buf, size_t n, size_t block_size, BitStream& bs){
    // 基准情形
    if (n <= 1) return;

    // 预计算步长（字节单位）
    const size_t stride_bytes = n/2 * block_size;


    // 递归调用：
    applyCompactionRec_inner<oswap_style>(current_buf, n/2, block_size, bs);

    applyCompactionRec_inner<oswap_style>(current_buf + (n/2) * block_size, n/2, block_size, bs);


    unsigned char* ptr_i = current_buf;                 // 左半区游标
    unsigned char* ptr_j = current_buf + (n/2) * block_size;  // 右半区游标

    for (size_t k = 0; k < n/2; ++k)
    {
        bool c = bs.next();

        oswap_buffer<oswap_style>(ptr_i, ptr_j, block_size,c);

        ptr_i += block_size;
        ptr_j += block_size;
    }
}

template <OSwap_Style oswap_style>
inline void applyCompactionRec(unsigned char* buf, size_t n, size_t block_size, const std::vector<uint8_t>& ctrlBits){

    FOAV_SAFE2_CNTXT(Butterfly_RecApply, n, block_size)
    if (!buf || n < 2) return;

    // Padding 处理 (保持原逻辑)
    const uint32_t Npad = nextPow2_u32((uint32_t)n);
    if (Npad < 2u) return;

    // 计算期望的控制位大小
    uint32_t l = 0;
    for (uint32_t tmp = Npad; tmp > 1u; tmp >>= 1u) ++l;

    const size_t totalSwitches = (size_t)(Npad / 2) * (size_t)l;
    const size_t expectedSize = (totalSwitches + 7) / 8;
    
    if (ctrlBits.size() != expectedSize) {
        printf("applyCompactionRec: ctrlBits.size=%zu expected=%zu\n",
               ctrlBits.size(), expectedSize);
        return;
    }

    // 内存准备 (Work Buffer Setup)
    unsigned char* work = buf;
    unsigned char* tmp = nullptr;
    
    // 如果不是2的幂，需要分配临时空间并 Padding
    if ((uint32_t)n != Npad) {
        tmp = (unsigned char*)malloc((size_t)Npad * block_size);
        if (!tmp) {
            printf("applyCompactionRec: malloc failed (Npad=%u)\n", Npad);
            return;
        }
        std::memcpy(tmp, buf, n * block_size);
        // Padding 部分清零（可选，视具体需求而定）
        std::memset(tmp + n * block_size, 0, ((size_t)Npad - n) * block_size);
        work = tmp;
    }

    // 执行核心算法
    BitStream bs(ctrlBits);

    applyCompactionRec_inner<oswap_style>(work, Npad, block_size, bs);

    // 清理与拷回
    if (tmp) {
        std::memcpy(buf, tmp, n * block_size);
        free(tmp);
    }
}


//迭代实现递归
template <OSwap_Style oswap_style>
static void applyCompactionRecIter_inner(unsigned char* root_buf, size_t root_len, size_t block_size, BitStream& bs){
    struct Frame {
        unsigned char* base; // 当前子数组起始地址
        size_t len;          // 子数组长度
        uint8_t state;       // 0=准备走左; 1=准备走右; 2=做merge并退出
    };

    // 用 vector 当栈，比 std::stack 更快一些
    std::vector<Frame> st;
    st.reserve(64); // log2(2^20)=20，64 足够覆盖很深的递归

    st.push_back(Frame{root_buf, root_len, 0});

    while (!st.empty()) {
        Frame& f = st.back();

        if (f.len <= 1) {
            st.pop_back();
            continue;
        }

        const size_t half = f.len >> 1;
        const size_t stride_bytes = half * block_size;

        if (f.state == 0) {
            // 先处理左半
            f.state = 1;
            st.push_back(Frame{f.base, half, 0});
            continue;
        }

        if (f.state == 1) {
            // 再处理右半
            f.state = 2;
            st.push_back(Frame{f.base + stride_bytes, half, 0});
            continue;
        }

        // f.state == 2: 左右都完成 -> 执行 merge（消耗 ctrl bits）
        unsigned char* ptr_i = f.base;
        unsigned char* ptr_j = f.base + stride_bytes;

        for (size_t k = 0; k < half; ++k) {
            const uint8_t c = bs.next();               // 与递归完全一致的读取顺序
            oswap_buffer<oswap_style>(ptr_i, ptr_j, (uint32_t)block_size, c);
            ptr_i += block_size;
            ptr_j += block_size;
        }

        st.pop_back();
    }
}

template <OSwap_Style oswap_style>
inline void applyCompactionRecIter(unsigned char* buf, size_t n, size_t block_size, const std::vector<uint8_t>& ctrlBits){
    FOAV_SAFE2_CNTXT(Butterfly_IterDFSApply, n, block_size)
    if (!buf || n < 2) return;

    const uint32_t Npad = nextPow2_u32((uint32_t)n);
    if (Npad < 2u) return;

    uint32_t l = 0;
    for (uint32_t tmp = Npad; tmp > 1u; tmp >>= 1u) ++l;

    const size_t totalSwitches = (size_t)(Npad / 2) * (size_t)l;
    const size_t expectedSize  = (totalSwitches + 7) / 8;
    if (ctrlBits.size() != expectedSize) {
        printf("applyCompactionRecIter: ctrlBits.size=%zu expected=%zu\n",
               ctrlBits.size(), expectedSize);
        return;
    }

    unsigned char* work = buf;
    unsigned char* tmp  = nullptr;

    if ((uint32_t)n != Npad) {
        tmp = (unsigned char*)malloc((size_t)Npad * block_size);
        if (!tmp) {
            printf("applyCompactionRecIter: malloc failed (Npad=%u)\n", Npad);
            return;
        }
        std::memcpy(tmp, buf, n * block_size);
        std::memset(tmp + n * block_size, 0, ((size_t)Npad - n) * block_size);
        work = tmp;
    }

    BitStream bs(ctrlBits);

    // 关键：使用 DFS 显式栈，访问顺序=递归
    applyCompactionRecIter_inner<oswap_style>(work, (size_t)Npad, block_size, bs);

    if (tmp) {
        std::memcpy(buf, tmp, n * block_size);
        free(tmp);
    }
}


//合并版
template <OSwap_Style oswap_style>
static void RecTripletKernel(unsigned char* buffer, size_t start, size_t len, size_t element_size, size_t meta_offset, uint32_t bit) {
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


//OR版
template <OSwap_Style oswap_style>
static void applyCompactionOR_inner2power(unsigned char* buf, size_t N, size_t block_size, BitStream& bs){
    if (N <= 1) return;
    if (N == 2) {
        const uint8_t sw = bs.next();
        //const uint8_t sw = 1;
        oswap_buffer<oswap_style>(buf, buf + block_size, (uint32_t)block_size, sw);
        return;
    }

    const size_t half = N >> 1;
    const size_t stride_bytes = half * block_size;

    // 先左后右（与 offline 生成一致）
    applyCompactionOR_inner2power<oswap_style>(buf,                 half, block_size, bs);
    applyCompactionOR_inner2power<oswap_style>(buf + stride_bytes,  half, block_size, bs);

    // 本层 merge：消费 half 个控制位
    unsigned char* L = buf;
    unsigned char* R = buf + stride_bytes;
    for (size_t i = 0; i < half; ++i) {
        const uint8_t sw = bs.next();
        //const uint8_t sw = 1;
        oswap_buffer<oswap_style>(L, R, (uint32_t)block_size, sw);
        L += block_size;
        R += block_size;
    }
}

template <OSwap_Style oswap_style>
static void applyCompactionOR_inner(unsigned char* buf, size_t N, size_t block_size, BitStream& bs){
    if (N <= 1) return;
        if (N == 2) {
            const uint8_t sw = bs.next();
            //const uint8_t sw = 1;
            oswap_buffer<oswap_style>(buf, buf + block_size, (uint32_t)block_size, sw);
            return;
        }

        // n1 = pow2_lt(N)   (最大2幂 < N)
        // n2 = N - n1
        const size_t n1 = pow2_lt_size(N);
        const size_t n2 = N - n1;

        unsigned char* L_ptr = buf;
        unsigned char* R_ptr = buf + n2 * block_size;

        // 递归处理左边 n2（general）
        applyCompactionOR_inner<oswap_style>(L_ptr, n2, block_size, bs);

        // 递归处理右边 n1（2power）
        applyCompactionOR_inner2power<oswap_style>(R_ptr, n1, block_size, bs);

        // merge：对齐 TightCompact_inner：

        unsigned char* R_suffix = buf + n1 * block_size;

        for (size_t i = 0; i < n2; ++i) {
            const uint8_t sw = bs.next();
            //const uint8_t sw = 1;
            oswap_buffer<oswap_style>(L_ptr, R_suffix, (uint32_t)block_size, sw);
            L_ptr += block_size;
            R_suffix += block_size;
        }
}

template <OSwap_Style oswap_style>
inline void applyCompactionOR(unsigned char* buf, size_t n, size_t block_size, const std::vector<uint8_t>& ctrlBits){
    FOAV_SAFE2_CNTXT(OR_applyCompaction, n, block_size)
    if (!buf || n < 2) return;

    BitStream bs(ctrlBits);

    applyCompactionOR_inner<oswap_style>(buf, n, block_size, bs);

}



template <OSwap_Style oswap_style>
static inline void OptApplyCompactionOR_inner2power(unsigned char* base, uint32_t len, size_t block_size, const uint8_t* __restrict c, size_t& t){
    if (len <= 1) return;

    if (len == 2) {
        uint8_t sw = c[t++];
        oswap_buffer<oswap_style>(base, base + block_size, (uint32_t)block_size, sw);
        return;
    }

    uint32_t half = len >> 1;
    size_t stride_bytes = (size_t)half * block_size;

    
    OptApplyCompactionOR_inner2power<oswap_style>(base, half, block_size, c, t);
    OptApplyCompactionOR_inner2power<oswap_style>(base + stride_bytes, half, block_size, c, t);

    // merge half swaps
    unsigned char* L = base;
    unsigned char* R = base + stride_bytes;
    for (uint32_t i = 0; i < half; ++i) {
        uint8_t sw = c[t++];
        oswap_buffer<oswap_style>(L, R, (uint32_t)block_size, sw);
        L += block_size;
        R += block_size;
    }
}

template <OSwap_Style oswap_style>
static inline void OptApplyCompactionOR_inner( unsigned char* base, uint32_t len, size_t block_size, const uint8_t* __restrict c, size_t& t){
    if (len <= 1) return;

    if (len == 2) {
        uint8_t sw = c[t++];
        oswap_buffer<oswap_style>(base, base + block_size, (uint32_t)block_size, sw);
        return;
    }

    // IMPORTANT: if len is pow2, do pure 2pow apply (avoids n2=0 general semantics)
    if (is_pow2_u32(len)) {
        OptApplyCompactionOR_inner2power<oswap_style>(base, len, block_size, c, t);
        return;
    }

    uint32_t n1 = pow2_le_u32(len);   // largest pow2 <= len
    uint32_t n2 = len - n1;

    OptApplyCompactionOR_inner<oswap_style>(base, n2, block_size, c, t);

    unsigned char* Rptr = base + (size_t)n2 * block_size;
    OptApplyCompactionOR_inner2power<oswap_style>(Rptr, n1, block_size, c, t);

    unsigned char* L = base;
    unsigned char* R_suffix = base + (size_t)n1 * block_size;

    for (uint32_t i = 0; i < n2; ++i) {
        uint8_t sw = c[t++];
        oswap_buffer<oswap_style>(L, R_suffix, (uint32_t)block_size, sw);
        L += block_size;
        R_suffix += block_size;
    }
}

/*
template <OSwap_Style oswap_style>
inline void OptApplyCompactionOR( unsigned char* buf, size_t n, size_t block_size, const std::vector<uint8_t>& ctrlBits){
    FOAV_SAFE2_CNTXT(OR_applyCompaction, n, block_size)
    if (!buf || n < 2) return;

    uint32_t N = (uint32_t)n;
    if (N < 2u) return;

    size_t totalSwitches = count_or_switches(N);
    if (ctrlBits.size() != totalSwitches) {
        printf("applyCompactionOR: ctrlBits.size=%zu expected=%zu\n",
               ctrlBits.size(), totalSwitches);
        return;
    }

    const uint8_t* __restrict c = ctrlBits.data();
    size_t t = 0;

    // Root: if N is pow2, go 2pow; else go general
    if (is_pow2_u32(N)) {
        OptApplyCompactionOR_inner2power<oswap_style>(buf, N, block_size, c, t);
    } else {
        OptApplyCompactionOR_inner<oswap_style>(buf, N, block_size, c, t);
    }

    // Debug: ensure all bits consumed (optional)
    // assert(t == totalSwitches);
}
*/

template <OSwap_Style oswap_style>
inline void OptApplyCompactionOR(unsigned char* buf, size_t n, size_t block_size, const std::vector<uint8_t>& ctrlBits){
    FOAV_SAFE2_CNTXT(OR_applyCompaction, n, block_size)
    if (!buf || n < 2) return;

    uint32_t N = (uint32_t)n;
    if (N < 2u) return;

    size_t totalSwitches = count_or_switches(N);
    if (ctrlBits.size() != totalSwitches) {
        printf("applyCompactionOR: ctrlBits.size=%zu expected=%zu\n",
               ctrlBits.size(), totalSwitches);
        return;
    }

    const uint8_t* __restrict c = ctrlBits.data();
    size_t t = 0;

    enum Kind : uint8_t { GEN = 0, POW2 = 1 };

    struct Frame {
        unsigned char* base;
        uint32_t len;
        uint8_t state; // 0,1,2
        Kind kind;
        // 仅 GEN 用：
        uint32_t n1;
        uint32_t n2;
    };

    std::vector<Frame> st;
    st.reserve(64);

    // root kind：若本身是 pow2，直接 POW2，否则 GEN
    Kind rootKind = is_pow2_u32(N) ? POW2 : GEN;
    st.push_back(Frame{buf, N, 0u, rootKind, 0u, 0u});

    while (!st.empty()) {
        Frame& f = st.back();

        if (f.len <= 1) { st.pop_back(); continue; }

        // len==2 统一处理（无论 GEN/POW2），消费 1 个 ctrl byte
        if (f.len == 2) {
            uint8_t sw = c[t++];
            oswap_buffer<oswap_style>(f.base, f.base + block_size, (uint32_t)block_size, sw);
            st.pop_back();
            continue;
        }

        // 若 kind==GEN 但 len 是 pow2，直接按 POW2 执行（避免 n2=0 offset 语义问题）
        if (f.kind == GEN && is_pow2_u32(f.len)) {
            f.kind = POW2;
            f.state = 0;
            continue;
        }

        if (f.kind == POW2) {
            uint32_t half = f.len >> 1;
            size_t stride_bytes = (size_t)half * block_size;

            if (f.state == 0) {
                f.state = 1;
                st.push_back(Frame{f.base, half, 0u, POW2, 0u, 0u});
                continue;
            }
            if (f.state == 1) {
                f.state = 2;
                st.push_back(Frame{f.base + stride_bytes, half, 0u, POW2, 0u, 0u});
                continue;
            }

            // merge half swaps
            unsigned char* L = f.base;
            unsigned char* R = f.base + stride_bytes;
            for (uint32_t i = 0; i < half; ++i) {
                uint8_t sw = c[t++];
                oswap_buffer<oswap_style>(L, R, (uint32_t)block_size, sw);
                L += block_size;
                R += block_size;
            }
            st.pop_back();
            continue;
        }

        // f.kind == GEN, len is NOT pow2 here
        if (f.state == 0) {
            uint32_t n1 = pow2_le_u32(f.len);
            uint32_t n2 = f.len - n1;
            f.n1 = n1;
            f.n2 = n2;
            f.state = 1;

            // prefix general
            st.push_back(Frame{f.base, n2, 0u, GEN, 0u, 0u});
            continue;
        }

        if (f.state == 1) {
            f.state = 2;
            // suffix pow2
            unsigned char* Rptr = f.base + (size_t)f.n2 * block_size;
            st.push_back(Frame{Rptr, f.n1, 0u, POW2, 0u, 0u});
            continue;
        }

        // state == 2: final merge (n2 swaps): L prefix vs suffix of R (start at base + n1*block_size)
        unsigned char* L = f.base;
        unsigned char* R_suffix = f.base + (size_t)f.n1 * block_size;
        for (uint32_t i = 0; i < f.n2; ++i) {
            uint8_t sw = c[t++];
            oswap_buffer<oswap_style>(L, R_suffix, (uint32_t)block_size, sw);
            L += block_size;
            R_suffix += block_size;
        }
        st.pop_back();
    }

    // 调试期可 assert(t == totalSwitches)
}

#endif // __BUTTERFLY_COMPACTION_TCC__