#include <stdio.h>
#include <stdalign.h>

static void DCT(const int * dataptr, int * outptr ) {
   int tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;
   tmp7 = dataptr[0] - dataptr[7];
   tmp0 = dataptr[0] - (tmp7 / 2);

   tmp6 = dataptr[1] - dataptr[6];
   tmp1 = dataptr[1] - (tmp6 / 2);

   tmp5 = dataptr[2] - dataptr[5];
   tmp2 = dataptr[2] - (tmp5 / 2);

   tmp4 = dataptr[3] - dataptr[4];
   tmp3 = dataptr[3] - (tmp4 / 2);
   //0-3
   tmp0 += tmp3;
   tmp3 = (tmp0 / 2)-tmp3;
   tmp1 += tmp2;
   tmp2 = (tmp1 / 2)-tmp2;

   outptr[0*8] = tmp0 = tmp0 + tmp1;
   outptr[4*8] = (tmp0 / 2)-tmp1;
   outptr[6*8] = tmp2 = tmp2-(tmp3*3/8);
   outptr[2*8] = (tmp2*3/8) + tmp3;
   //4-7
   tmp6 += (tmp5*3/8);
   tmp5 = (tmp6*5/8)-tmp5;

   tmp4 += tmp5;
   tmp5 = (tmp4/2)-tmp5;

   outptr[1*8] = tmp7 = tmp7+tmp6;
   tmp6 = (tmp7/2) - tmp6;
   outptr[7*8] = tmp4-(tmp7/8);

   outptr[5*8] = tmp5 = tmp5 + (tmp6*7/8);
   outptr[3*8] = tmp6-(tmp5/2);
}

static void IDCT(const int * dataptr, int * outptr) {
   int tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;
   tmp4 =  (dataptr[0] / 2) - dataptr[4];
   tmp0 =  dataptr[0] - tmp4;

   tmp2 = dataptr[2] - (dataptr[6]*3/8);
   tmp6 = dataptr[6] + (tmp2*3/8);

   tmp3 += (dataptr[5] / 2);
   tmp5 = dataptr[5] - (tmp3*7/8);

   tmp1 = dataptr[1];
   tmp7 = dataptr[7];

   //0-4-6-2print_matrix_type(matrix,"%d",8,1)
   tmp6 = (tmp4 / 2) - tmp6;
   tmp4 = tmp4 - tmp6;

   tmp2 = (tmp0 / 2) - tmp2;
   tmp0 = tmp0 - tmp2;

   //7-5-3-1
   tmp7 = tmp7 + (tmp1/8);

   tmp5 = (tmp7 / 2) - tmp5;
   tmp7 = tmp7 - tmp5;

   tmp3 = (tmp1 / 2) - tmp3;
   tmp1 = tmp1 - tmp3;

   tmp5 = (tmp3*5/8) - tmp5;
   tmp3 = tmp3 - (tmp5*3/8);
   //join
   outptr[3*8] = tmp2 = tmp2 + (tmp7 / 2);
   outptr[4*8] = tmp2 - tmp7;

   outptr[2*8] = tmp6 = tmp6 + (tmp5 / 2);
   outptr[5*8] = tmp6 - tmp5;

   outptr[1*8] = tmp4 = tmp4 + (tmp3 / 2);
   outptr[6*8] = tmp4 - tmp3;

   outptr[0*8] = tmp0 = tmp0 + (tmp1 / 2);
   outptr[7*8] = tmp0 - tmp1;
}



void dct8x8(const int*in, int*out) {
   alignas(32) int tmp[64];
   for (int i=0; i<8; i++) {
      DCT(&in[i*8],&tmp[i]);
   }
   for (int i=0; i<8; i++) {
      DCT(&tmp[i*8],&out[i]);
   }
}

void idct8x8(const int*in, int*out) {
   alignas(32) int tmp[64];
   for (int i=0; i<8; i++) {
      IDCT(&in[i],&tmp[i*8]);
   }
   for (int i=0; i<8; i++) {
      IDCT(&tmp[i],&out[i*8]);
   }
}
