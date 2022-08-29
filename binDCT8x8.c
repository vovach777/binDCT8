#include <stdalign.h>
static void _row8_bdct(int dataptr[8], int dataptr_[8*8] ){
   int tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp10,tmp11,tmp12, tmp13;
   tmp0 = dataptr[0] + dataptr[7];
   tmp7 = dataptr[0] - dataptr[7];
   tmp1 = dataptr[1] + dataptr[6];
   tmp6 = dataptr[1] - dataptr[6];
   tmp2 = dataptr[2] + dataptr[5];
   tmp5 = dataptr[2] - dataptr[5];
   tmp3 = dataptr[3] + dataptr[4];
   tmp4 = dataptr[3] - dataptr[4];
   tmp5 = tmp5 - (tmp6 >> 1);
   tmp6 = tmp5 - (tmp5>>2) + tmp6;
   tmp5 = (tmp6>>1)   - tmp5;
   tmp10 = tmp0 + tmp3;
   tmp13 = tmp0 - tmp3;
   tmp11 = tmp1 + tmp2;
   tmp12 = tmp1 - tmp2;
   dataptr_[0*8] = tmp10 = tmp10+tmp11;
   dataptr_[4*8] = ((tmp10>>1) - tmp11)<<1;
   dataptr_[6*8] = tmp12 = (tmp13>>1)-tmp12;
   dataptr_[2*8] = (tmp13 - (tmp12>>1))<<1;
   tmp10 = tmp4+tmp5;
   tmp11 = tmp4-tmp5;
   tmp12 = tmp7 - tmp6;
   tmp13 = tmp7 + tmp6;
   dataptr_[7*8] = tmp10 = ((tmp13>>2)-tmp10)>>1;
   dataptr_[1*8] = (tmp13 - (tmp10>>2))<<1; //scale x2
   dataptr_[5*8] = tmp11 = tmp12 + tmp11;
   dataptr_[3*8] = (tmp12 - (tmp11>>1))<<1; //scale x2
}


static void dctrows(int in_data[8*8], int out_data[8*8]) {
  for (int i=0;i<8;i++) {
    _row8_bdct(&in_data[i*8],  &out_data[i]);
  }
}

void bdct_8x8(int data[8*8]) {
alignas(32) int tmp[8*8];
   dctrows(data,tmp);
   dctrows(tmp,data);
}
