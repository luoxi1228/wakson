#ifndef __BUTTERFLYNETWORK_HPP__
#define __BUTTERFLYNETWORK_HPP__

#include <cstdint>
#include <vector>
#include <cstring>

#include "../oasm_lib.h"
#include "../ObliviousPrimitives.hpp"
#include "../Enclave_globals.h"
#include "../utils.hpp"

/*
  ButterflyNetwork: canonical logN-stage butterfly of 2x2 switches.

  Offline:
    setControlBits(sel, Npad):
      - sel length must be exactly Npad (Npad is power of two)
      - computes and stores control bit for every 2x2 switch at every stage
      - applies the same OSwap semantics on internal packets to propagate tags

  Online:
    apply<oswap_style>(buf, block_size):
      - uses stored ctrlBits_ to OSwap data items in-place, stage by stage
*/

class ButterflyNetwork {
public:
  explicit ButterflyNetwork(uint32_t Npow2);

  // Offline: sel must be of length Ntotal_ and each entry in {0,1}
  void setControlBits(const uint8_t* sel);

  // Online
  template <OSwap_Style oswap_style>
  inline void apply(unsigned char* buf, size_t block_size) const;

  // Debug / info
  uint32_t getSize() const { return Ntotal_; }  // 网络总大小
  uint32_t getNumStages() const { return numStages_; } // 网络级数 = log2(N)
  uint32_t getSwitchesPerStage() const { return switchesPerStage_; } // 每级开关数 = N/2

  // 判断是否为 2 的幂次
  static inline bool isPowerOfTwo(uint32_t x) { 
    return x && ((x & (x - 1)) == 0); 
  }

  // 计算大于等于 x 的最小 2 的幂次
  static inline uint32_t nextPow2(uint32_t x) {
    if (x <= 1) return 1;
    // round up to next power of two
    --x;
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    return x + 1;
  }

private:
  uint32_t Ntotal_;            // must be power-of-two
  uint32_t numStages_;         // log2(Ntotal_)
  uint32_t switchesPerStage_;  // N/2

  // stage-major: ctrlBits_[stage * switchesPerStage_ + sw] ∈ {0,1}
  std::vector<uint8_t> ctrlBits_;

  inline uint8_t getCtrl_(uint32_t stage, uint32_t sw) const {
    return ctrlBits_[(size_t)stage * (size_t)switchesPerStage_ + (size_t)sw];
  }
  inline void setCtrl_(uint32_t stage, uint32_t sw, uint8_t b) {
    ctrlBits_[(size_t)stage * (size_t)switchesPerStage_ + (size_t)sw] = (uint8_t)(b & 1u);
  }

  // Online: apply one stage
  template <OSwap_Style oswap_style>
  inline void applyStage_(unsigned char* buf, size_t block_size, uint32_t stage) const;
};

// ------------------------
// Online apply (template)
// ------------------------
template <OSwap_Style oswap_style>
inline void ButterflyNetwork::apply(unsigned char* buf, size_t block_size) const {
  //安全检查
  FOAV_SAFE_CNTXT(Butterfly_AP, Ntotal_)
  if (Ntotal_ < 2) return;   //少于两个元素不需要交换

  //逐层应用交换
  for (uint32_t s = 0; s < numStages_; ++s) {
    applyStage_<oswap_style>(buf, block_size, s);
  }
}

template <OSwap_Style oswap_style>
inline void ButterflyNetwork::applyStage_(unsigned char* buf,
                                         size_t block_size,
                                         uint32_t stage) const {
  //安全检查
  FOAV_SAFE_CNTXT(Butterfly_APStage, Ntotal_)
  FOAV_SAFE_CNTXT(Butterfly_APStage, stage)

  //计算当前网络层级的参数
  const uint32_t d = stage;      //层级编号
  const uint32_t stride = (1u << d);  //跨度 = 2^d

  // 遍历所有交换对
  for (uint32_t b = 0; b < Ntotal_; b += (stride << 1)) {
    for (uint32_t k = 0; k < stride; ++k) {
      //计算交换索引
      const uint32_t i = b + k;
      const uint32_t j = b + k + stride;

      //计算开关的索引
      const uint32_t sw = (b >> 1) + k; // 0..N/2-1 stage-local switch index
      const uint8_t c = getCtrl_(stage, sw);  //获取控制位

      oswap_buffer<oswap_style>(
          buf + (size_t)block_size * (size_t)i,
          buf + (size_t)block_size * (size_t)j,
          (uint32_t)block_size,
          c);
    }
  }
}


void OblivButterflyCompact(unsigned char* buffer, uint32_t n, size_t block_size, enc_ret* ret);

void DecryptAndOblivButterflyCompact(unsigned char* encrypted_buffer, uint32_t n,
                                     size_t encrypted_block_size,
                                     unsigned char* result_buffer,
                                     enc_ret* ret);

#endif // __BUTTERFLYNETWORK_HPP__
