#pragma once

typedef char BitVectorType;

class BitVector
{
public :
  BitVector(int n) {
    int n_ = (n + sizeof(BitVectorType) - 1)/sizeof(BitVectorType);
    bv_ = new BitVectorType[n_];
#pragma omp parallel for
    for (int i = 0; i < n_; ++i) {
      bv_[i] = 0;
    }
  }

  ~BitVector() { delete[] bv_; }

  void set(int i) {
    bv_[getIndexOf_(i)] |= getMaskOf_(i);
  }

  bool get(int i) const {
    return bv_[getIndexOf_(i)] & getMaskOf_(i);
  }

  bool testAndSet(int i) {
    if (!get(i)) {
      BitVectorType mask = getMaskOf_(i);
      BitVectorType prev = __sync_fetch_and_or(bv_ + getIndexOf_(i), mask);
      return !(prev & mask);
    }
    else {
      return false;
    }
  }

  bool atomicClear(int i) {
    __sync_fetch_and_and(bv_ + getIndexOf_(i), ~getMaskOf_(i));
  }

private :
  typedef char BitVectorType;

  static int getIndexOf_(int i) { return i/sizeof(BitVectorType); }
  static int getBitIndexOf_(int i) { return i%sizeof(BitVectorType); }
  static BitVectorType getMaskOf_(int i) { return 1 << getBitIndexOf_(i); }

  BitVectorType *bv_;
};
