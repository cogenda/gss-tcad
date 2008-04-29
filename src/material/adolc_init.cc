#include "adolc.h"
unsigned int adtl::AutoDScalar::numdir = 9;

extern "C"
{
  void  set_ad_number(const unsigned int p)
  {
    adtl::AutoDScalar::numdir=p;
  }
}
