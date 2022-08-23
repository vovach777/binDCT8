void binDCT8_transform(int in_out[8])
{
    int tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, tmp17;

        //stage 1
        tmp0 = in_out[0] + in_out[7];
        tmp7 = in_out[0] - in_out[7];
        tmp1 = in_out[1] + in_out[6]; 
        tmp6 = in_out[1] - in_out[6];
        tmp2 = in_out[2] + in_out[5];
        tmp5 = in_out[2] - in_out[5];
        tmp3 = in_out[3] + in_out[4];
        tmp4 = in_out[3] - in_out[4];

        //stage 2
        //tmp16 = ((tmp5*3)>>3) + tmp6;
        //tmp15 = ((tmp16*5)>>3) - tmp5;
        tmp16 = ((tmp5<<2)-tmp5>>3) + tmp6;
        tmp15 = ((tmp16<<2)+tmp16>>3) - tmp5;

        //stage 3
        tmp10 = tmp0 + tmp3;
        tmp13 = tmp0 - tmp3;
        tmp11 = tmp1 + tmp2;
        tmp12 = tmp1 - tmp2;

        tmp14 = tmp4 + tmp15;
        tmp15 = tmp4 - tmp15;


        tmp0 = tmp16;
        tmp16 = tmp7 - tmp16;
        tmp17 = tmp0 + tmp7;

        //stage 4
        tmp14 = (tmp17 >> 3) - tmp14;

        tmp10 = tmp10 + tmp11;
        tmp11 = (tmp10 >> 1) - tmp11;

        //tmp12 = ((tmp13*3)>>3) - tmp12;
        //tmp13 = ((tmp12*3)>>3) + tmp13;
        tmp12 = ((tmp13 << 2)-tmp13>>3) - tmp12;
        tmp13 = ((tmp12 << 2)-tmp12>>3) + tmp13;

        //tmp15 = ((tmp16*7)>>3) + tmp15;
        tmp15 = ((tmp16 << 3)-tmp16>>3) + tmp15;
        tmp16 = (tmp15>>1) - tmp16;

        //stage 5
        in_out[0] = tmp10;
        in_out[4] = tmp11;
        in_out[6] = tmp12;
        in_out[2] = tmp13;
        in_out[7] = tmp14;
        in_out[5] = tmp15;
        in_out[3] = tmp16;
        in_out[1] = tmp17;
   
}
