#include "ButterflyNetwork.hpp"

#include <cstdlib>   // malloc/free
#include <cstdio>    // printf


// Constructor: ButterflyNetwork(Npow2)
ButterflyNetwork::ButterflyNetwork(uint32_t Npow2)
  : Ntotal_(Npow2),
    numStages_(0),
    switchesPerStage_(0),
    ctrlBits_() {

  FOAV_SAFE_CNTXT(Butterfly_CTOR, Ntotal_)

  // Enclave-friendly: no throw
  if (!isPowerOfTwo(Ntotal_) || Ntotal_ < 1) {
    printf("ButterflyNetwork: N must be power-of-two and >= 1 (got %u)\n", Ntotal_);
    // Put object into a safe empty state
    Ntotal_ = 1;
    numStages_ = 0;
    switchesPerStage_ = 0;
    ctrlBits_.clear();
    return;
  }

  uint32_t tmp = Ntotal_;
  while (tmp > 1) {
    tmp >>= 1;
    ++numStages_;
  }

  switchesPerStage_ = Ntotal_ / 2;
  ctrlBits_.assign((size_t)numStages_ * (size_t)switchesPerStage_, 0);
}

// Internal packet for offline control-bit propagation.
// Use fixed 8-byte layout so OSWAP_8 is always valid.
struct BNPacket {
  uint32_t value; // 0/1
  uint32_t tag;   // rank tag
};
static_assert(sizeof(BNPacket) == 8, "BNPacket must be exactly 8 bytes");

// Offline: compute control bits from sel[0..Ntotal_-1]
void ButterflyNetwork::setControlBits(const uint8_t* sel) {
  FOAV_SAFE_CNTXT(Butterfly_SetCB, Ntotal_)

  if (Ntotal_ < 2 || numStages_ == 0) return;

  // 1) Init packets P[i] = <value, tag>
  std::vector<BNPacket> P((size_t)Ntotal_);
  for (uint32_t i = 0; i < Ntotal_; ++i) {
    P[(size_t)i].value = (uint32_t)(sel[i] & 1u);
    P[(size_t)i].tag   = 0;
  }

  // 2) Ranking: tag = rank * value; rank += value
  uint32_t rank = 0;
  for (uint32_t i = 0; i < Ntotal_; ++i) {
    const uint32_t v = P[(size_t)i].value;
    P[(size_t)i].tag = rank * v;
    rank += v;
  }

  // 3) Stage-by-stage: compute c and OSwap packets
  for (uint32_t d = 0; d < numStages_; ++d) {
    const uint32_t stride = (1u << d);

    for (uint32_t b = 0; b < Ntotal_; b += (stride << 1)) {
      for (uint32_t k = 0; k < stride; ++k) {
        const uint32_t i = b + k;
        const uint32_t j = b + k + stride;

        const uint32_t bit_i = (P[(size_t)i].tag >> d) & 1u;
        const uint32_t bit_j = (P[(size_t)j].tag >> d) & 1u;

        const uint32_t down_i = P[(size_t)i].value & bit_i;
        const uint32_t up_j   = P[(size_t)j].value & (bit_j ^ 1u);
        const uint8_t  c      = (uint8_t)((down_i | up_j) & 1u);

        const uint32_t sw = (b >> 1) + k;
        setCtrl_(d, sw, c);

        // Propagate tags/values through same switch pattern
        oswap_buffer<OSWAP_8>(
            (unsigned char*)&P[(size_t)i],
            (unsigned char*)&P[(size_t)j],
            8u,
            c);
      }
    }
  }
}

// -----------------------------
// Wrapper: OblivButterflyCompact
// -----------------------------
void OblivButterflyCompact(unsigned char* buffer, uint32_t n, size_t block_size, enc_ret* ret) {
  FOAV_SAFE_CNTXT(OblivButterflyCompact, n)
  if (n < 2) return;

  const uint32_t Npad = ButterflyNetwork::nextPow2(n);

  long t1, t2;

  // Build sel padded to Npad
  std::vector<uint8_t> sel((size_t)Npad, 0);
  // for (uint32_t i = 0; i < n; ++i) {
  //   sel[(size_t)i] = (*(uint8_t*)(buffer + (size_t)block_size * (size_t)i)) & 1u;
  // }
  for (uint32_t i = 0; i < n; ++i) {
    uint8_t r = 0;
    getRandomBytes(&r, 1);     
    sel[(size_t)i] = (uint8_t)(r & 1u);
  }

  // pad 的部分保持 0（虚拟元素不选中）
  for (uint32_t i = n; i < Npad; ++i) {
    sel[(size_t)i] = 0;
  }

  // If already power-of-two, we can do in-place
  if (Npad == n) {
    ocall_clock(&t1);
    ButterflyNetwork bnet(Npad);
    bnet.setControlBits(sel.data());
    ocall_clock(&t2);
    if (ret) ret->control_bits_time = ((double)(t2 - t1)) / 1000.0;

#ifdef COUNT_OSWAPS
    if (ret) ret->OSWAP_cb = OSWAP_COUNTER;
    OSWAP_COUNTER = 0;
#endif

    ocall_clock(&t1);
    if (block_size == 4) {
      bnet.apply<OSWAP_4>(buffer, block_size);
    } else if (block_size == 8) {
      bnet.apply<OSWAP_8>(buffer, block_size);
    } else if (block_size == 12) {
      bnet.apply<OSWAP_12>(buffer, block_size);
    } else if (block_size % 16 == 0) {
      bnet.apply<OSWAP_16X>(buffer, block_size);
    } else {
      bnet.apply<OSWAP_8_16X>(buffer, block_size);
    }
    ocall_clock(&t2);
    if (ret) ret->apply_perm_time = ((double)(t2 - t1)) / 1000.0;

#ifdef COUNT_OSWAPS
    if (ret) ret->OSWAP_ap = OSWAP_COUNTER;
#endif
    return;
  }

  // Otherwise: allocate temp buffer of Npad records, copy, pad with zeros
  unsigned char* tmp = (unsigned char*)malloc((size_t)Npad * block_size);
  if (!tmp) {
    printf("OblivButterflyCompact: malloc failed for tmp buffer\n");
    return;
  }
  // copy real records
  std::memcpy(tmp, buffer, (size_t)n * block_size);
  // pad dummy records (all-zero)
  std::memset(tmp + (size_t)n * block_size, 0, ((size_t)Npad - (size_t)n) * block_size);

  // Offline
  ocall_clock(&t1);
  ButterflyNetwork bnet(Npad);
  bnet.setControlBits(sel.data());
  ocall_clock(&t2);
  if (ret) ret->control_bits_time = ((double)(t2 - t1)) / 1000.0;

#ifdef COUNT_OSWAPS
  if (ret) ret->OSWAP_cb = OSWAP_COUNTER;
  OSWAP_COUNTER = 0;
#endif

  // Online
  ocall_clock(&t1);
  if (block_size == 4) {
    bnet.apply<OSWAP_4>(tmp, block_size);
  } else if (block_size == 8) {
    bnet.apply<OSWAP_8>(tmp, block_size);
  } else if (block_size == 12) {
    bnet.apply<OSWAP_12>(tmp, block_size);
  } else if (block_size % 16 == 0) {
    bnet.apply<OSWAP_16X>(tmp, block_size);
  } else {
    bnet.apply<OSWAP_8_16X>(tmp, block_size);
  }
  ocall_clock(&t2);
  if (ret) ret->apply_perm_time = ((double)(t2 - t1)) / 1000.0;

#ifdef COUNT_OSWAPS
  if (ret) ret->OSWAP_ap = OSWAP_COUNTER;
#endif

  // Copy back only first n records (same external footprint)
  std::memcpy(buffer, tmp, (size_t)n * block_size);
  free(tmp);
}

void DecryptAndOblivButterflyCompact(unsigned char* encrypted_buffer,
                                    uint32_t n,
                                    size_t encrypted_block_size,
                                    unsigned char* result_buffer,
                                    enc_ret* ret) {
  long t1, t2;

  unsigned char* decrypted_buffer = NULL;
  size_t decrypted_block_size =
      decryptBuffer(encrypted_buffer, (uint64_t)n, encrypted_block_size, &decrypted_buffer);

  ocall_clock(&t1);
  PRB_pool_init(1);

  OblivButterflyCompact(decrypted_buffer, n, decrypted_block_size, ret);

  ocall_clock(&t2);
  if (ret) ret->ptime = ((double)(t2 - t1)) / 1000.0;

#ifdef COUNT_OSWAPS
  if (ret) ret->OSWAP_count = OSWAP_COUNTER;
#endif

  encryptBuffer(decrypted_buffer, (uint64_t)n, decrypted_block_size, result_buffer);

  PRB_pool_shutdown();
  free(decrypted_buffer);
}
