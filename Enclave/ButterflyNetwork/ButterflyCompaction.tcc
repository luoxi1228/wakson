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
                            const std::vector<uint8_t>& ctrlBits) {
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

  const size_t expected = ((size_t)Npad / 2) * (size_t)l;
  if (ctrlBits.size() != expected) {
    printf("applyCompaction: ctrlBits.size=%zu expected=%zu (n=%zu Npad=%u l=%u)\n",
           ctrlBits.size(), expected, n, Npad, l);
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

        const uint8_t c = ctrlBits[t++];

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

#endif // __BUTTERFLY_COMPACTION_TCC__
