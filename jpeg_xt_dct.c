#define preshift 0
#define T int
#define deadzone 0
#define optimize 0

/// Multiplications by constants
#define FRACT_BITS 12
#define ROUND(x) (((x) + ((1 << FRACT_BITS) >> 1)) >> FRACT_BITS)
                                                            //                    109876543210
                                                            // Multiply by 403  = 000110010011
#define pmul_tan1(x) (t = (x) + ((x) << 1),t = t + ((x) << 4) + (t << 7),ROUND(t))

                                                            // Multiply by 1243 = 010011011011
#define pmul_tan3(x) (t = (x) + ((x) << 1),t = t + (t << 3) + (t << 6) + ((x) << 10),ROUND(t))

                                                            // Multiply by 1697 = 011010100001
#define pmul_tan4(x) (t = (x) + ((x) << 5) + ((x) << 7) + ((x) << 9) + ((x) << 10),ROUND(t))

                                                            // Multiply by 815  = 001100101111
#define pmul_tan2(x) (t = ((x) << 6) - ((x) << 4) - (x) + ((x) << 8) + ((x) << 9),ROUND(t))

                                                            // Multiply by 799  = 001100011111
#define pmul_sin1(x) (t = ((x) << 5) - (x) + ((x) << 8) + ((x) << 9),ROUND(t))

                                                            // Multiply by 2276 = 100011100100
#define pmul_sin3(x) (t = ((x) << 8) - ((x) << 5) + ((x) << 2) + ((x) << 11),ROUND(t))

                                                            // Multiply by 1567 = 011000011111
#define pmul_sin2(x) (t = ((x) << 5) - (x) + ((x) << 9) + ((x) << 10),ROUND(t))

                                                            // Multiply by 2896 = 101101010000
#define pmul_sin4(x) (t = (x) + ((x) << 2),t = ((x) << 4) + (t << 6) + (t << 9),ROUND(t))
///


void LiftingDCT_TransformBlock(const T *source,T *target)
{ 
  T *dpend,*dp;
  int band = 0;
  T t;
  //
  // Adjust the DC offset to the number of fractional bits.
  // The lifting DCT here has a scaling factor of 2 per dimension, and prescales data by FIX_BITS
  // three additional bits because we still need to divide by 8.
  //
  // Pass over columns
  dpend = target + 8;
  for(dp = target;dp < dpend;dp++,source++) {
    // Compute sqrt(2) T_8(0). This is the forward butterfly.
    T x0    = source[0 << 3] >> preshift;
    T x4    = source[7 << 3] >> preshift;
    x0     += pmul_tan4(x4);
    x4     -= pmul_sin4(x0);
    x0     += pmul_tan4(x4);
    x4      = -x4;
    T x1    = source[1 << 3] >> preshift;
    T x5    = source[6 << 3] >> preshift;
    x1     += pmul_tan4(x5);
    x5     -= pmul_sin4(x1);
    x1     += pmul_tan4(x5);
    x5      = -x5;
    T x2    = source[2 << 3] >> preshift;
    T x6    = source[5 << 3] >> preshift;
    x2     += pmul_tan4(x6);
    x6     -= pmul_sin4(x2);
    x2     += pmul_tan4(x6);
    x6      = -x6;
    T x3    = source[3 << 3] >> preshift;
    T x7    = source[4 << 3] >> preshift;
    x3     += pmul_tan4(x7);
    x7     -= pmul_sin4(x3);
    x3     += pmul_tan4(x7);
    x7      = -x7;
    // Compute the bold-Z vector from x0...x3 by T_4(0)
    T zb0    = x0 + pmul_tan4(x3);
    T zb2    = x3 - pmul_sin4(zb0);
    zb0     += pmul_tan4(zb2);
    zb2      = -zb2;
    T zb1    = x1 + pmul_tan4(x2);
    T zb3    = x2 - pmul_sin4(zb1);
    zb1     += pmul_tan4(zb3);
    zb3      = -zb3;
    // Compute z0,z1,z2 from the w-vector. This applies T_4(1) by two rotations, each of which consist of three shears.
    T z00    = pmul_tan1(x7)   + x4;
    T z01    = pmul_tan3(x6)   + x5;
    T z10    = -pmul_sin1(z00) + x7;
    T z11    = -pmul_sin3(z01) + x6;
    T z20    = pmul_tan1(z10)  + z00;
    T z21    = pmul_tan3(z11)  + z01;
    // Apply now T_8(0,1,0,0).
    // The output vector is now: zb0..zb3,z20,z21,z11,-z10.
    // This is called x^(2) in the bommer paper.
    // Compute again the bold Z vector as C_II oplus C_II. This is the lower half of the T_8(0,1,0,0) matrix.
    T zc0    = z20 + pmul_tan4(z21);
    T zc1    = z21 - pmul_sin4(zc0);
    zc0     += pmul_tan4(zc1);
    zc1      = -zc1;
    T zc3    = z11 + pmul_tan4(z10);
    T zc2    = z10 - pmul_sin4(zc3);
    zc3     += pmul_tan4(zc2);
    zc2      = -zc2;
    // Compute z0,z1,z2 from the w-vector. This is the upper half of the T_8(0,1,0,0) matrix and rotates by Pi/4 and Pi/8.
    z00 = pmul_tan4(zb1)  + zb0; // rotate zc0,zc1 by Pi/4
    z01 = pmul_tan2(zb3)  + zb2; // rotate zc2,zc3 by Pi/8
    z10 = -pmul_sin4(z00) + zb1;
    z11 = -pmul_sin2(z01) + zb3;
    z20 = pmul_tan4(z10)  + z00;
    z21 = pmul_tan2(z11)  + z01;
    // Compute the x3 vector: z20,-z10,z21,-z11,zc0,zc1,zc2,zc3
    // Not needed here, instead go ahead directly.
    // Instead, compute x4 directly. The first four comonets
    // are identical to x3, x47 = x63
    // This is the I_4 part of the last matrix.
    // Upper part of A_4(1) is identity. The matrix notation in the paper has here zc1,zc2 interchanged.
    T z0  = pmul_tan4(zc3) + zc1;
    T z1  = -pmul_sin4(z0) + zc3;
    T x45 = pmul_tan4(z1)  + z0; // rot by Pi/4
    // Output permutation by the B_8 matrix.
    dp[0 << 3] = z20;
    dp[1 << 3] = zc0;
    dp[2 << 3] = z21;
    dp[3 << 3] = -z1;
    dp[4 << 3] = -z10;
    dp[5 << 3] = x45;
    dp[6 << 3] = -z11;
    dp[7 << 3] = zc2;
  }
  //
  // Pass over rows and quantize, remove the dc shift.
  dpend = target + (8 << 3);
  for(dp = target;dp < dpend;dp += 8) {
    // Step 1: Compute sqrt(2) T_8(0). This is the forward butterfly.
    T x0    = dp[0];
    T x4    = dp[7];
    x0     += pmul_tan4(x4);
    x4     -= pmul_sin4(x0);
    x0     += pmul_tan4(x4);
    x4      = -x4;
    T x1    = dp[1];
    T x5    = dp[6];
    x1     += pmul_tan4(x5);
    x5     -= pmul_sin4(x1);
    x1     += pmul_tan4(x5);
    x5      = -x5;
    T x2    = dp[2];
    T x6    = dp[5];
    x2     += pmul_tan4(x6);
    x6     -= pmul_sin4(x2);
    x2     += pmul_tan4(x6);
    x6      = -x6;
    T x3    = dp[3];
    T x7    = dp[4];
    x3     += pmul_tan4(x7);
    x7     -= pmul_sin4(x3);
    x3     += pmul_tan4(x7);
    x7      = -x7;
    // Compute the bold-Z vector from x0...x3 by T_4(0) 
    T zb0    = x0 + pmul_tan4(x3);
    T zb2    = x3 - pmul_sin4(zb0);
    zb0     += pmul_tan4(zb2);
    zb2      = -zb2;
    T zb1    = x1 + pmul_tan4(x2);
    T zb3    = x2 - pmul_sin4(zb1);
    zb1     += pmul_tan4(zb3);
    zb3      = -zb3;
    // Compute z0,z1,z2 from the w-vector. This applies T_4(1) by two rotations, each of which consist of three shears.
    T z00    = pmul_tan1(x7)   + x4;
    T z01    = pmul_tan3(x6)   + x5;
    T z10    = -pmul_sin1(z00) + x7;
    T z11    = -pmul_sin3(z01) + x6;
    T z20    = pmul_tan1(z10)  + z00;
    T z21    = pmul_tan3(z11)  + z01;
    // The output vector is now: zb0..zb3,z20,z21,z11,-z10.
    // This is called x^(2) in the bommer paper.
    // Compute again the bold Z vector as C_II oplus C_II
    T zc0    = z20 + pmul_tan4(z21);
    T zc1    = z21 - pmul_sin4(zc0);
    zc0     += pmul_tan4(zc1);
    zc1      = -zc1;
    T zc3    = z11 + pmul_tan4(z10);
    T zc2    = z10 - pmul_sin4(zc3);
    zc3     += pmul_tan4(zc2);
    zc2      = -zc2;
    // Compute z0,z1,z2 from the w-vector
    z00 = pmul_tan4(zb1)  + zb0;
    z01 = pmul_tan2(zb3)  + zb2;
    z10 = -pmul_sin4(z00) + zb1;
    z11 = -pmul_sin2(z01) + zb3;
    z20 = pmul_tan4(z10)  + z00;
    z21 = pmul_tan2(z11)  + z01;
    // Compute the x3 vector: z20,-z10,z21,-z11,zc0,zc1,zc2,zc3
    // Not needed here, instead go ahead directly.
    // Instead, compute x4 directly. The first four componets
    // are identical to x3, zc2 = x63
    // Upper part of last matrix is identity
    // Upper part of A_4(1) is identity
    T z0  = pmul_tan4(zc3) + zc1;
    T z1  = -pmul_sin4(z0) + zc3;
    T x45 = pmul_tan4(z1)  + z0;
    // Output permutation.
    dp[0] = z20;
    dp[1] = zc0;
    dp[2] = z21;
    dp[3] = -z1;
    dp[4] = -z10;
    dp[5] = x45;
    dp[6] = -z11;
    dp[7] = zc2 ;

    band    += 8;
  }
}
///

/// LiftingDCT::InverseTransformBlock
// Run the inverse DCT on an 8x8 block reconstructing the data.

void LiftingDCT_InverseTransformBlock(const T *source, T *target)                                                                 
{
    T t;
    T *dp;
    T *dpend = target + (8 << 3);
    for(dp = target;dp < dpend;dp += 8, source += 8) {
      // Inverse output permutation.
      T z20 =  source[0];
      T zc0 =  source[1];
      T z21 =  source[2];
      T  z1 = -source[3];
      T z10 = -source[4];
      T x45 =  source[5];
      T z11 = -source[6];
      T zc2 =  source[7];
      // rotation by 45 degrees of x45,x46 to zc1,zc3.
      T z0  = x45 - pmul_tan4(z1);
      T zc3 = z1  + pmul_sin4(z0);
      T zc1 = z0  - pmul_tan4(zc3);
      // Next rotation pair.
      T z00 = z20 - pmul_tan4(z10);
      T z01 = z21 - pmul_tan2(z11);
      T zb1 = z10 + pmul_sin4(z00);
      T zb3 = z11 + pmul_sin2(z01);
      T zb0 = z00 - pmul_tan4(zb1);
      T zb2 = z01 - pmul_tan2(zb3);
      // Small butterfly. 
      zc1      = -zc1;
      zc0     -= pmul_tan4(zc1);
      z21      = zc1 + pmul_sin4(zc0);
      z20      = zc0 - pmul_tan4(z21);
      zc2      = -zc2;
      zc3     -= pmul_tan4(zc2);
      z10      = zc2 + pmul_sin4(zc3);
      z11      = zc3 - pmul_tan4(z10);
      // Rotation by 3Pi/16 and 1Pi/16.
      z00      = z20 - pmul_tan1(z10);
      z01      = z21 - pmul_tan3(z11);
      T x7  = z10 + pmul_sin1(z00);
      T x6  = z11 + pmul_sin3(z01);
      T x4  = z00 - pmul_tan1(x7);
      T x5  = z01 - pmul_tan3(x6);
      // Small butterfly again to compute x0,x1,x2,x3, again only inverted up to a scale factor of 2.
      zb2      = -zb2;
      zb0     -= pmul_tan4(zb2);
      T x3  = zb2 + pmul_sin4(zb0);
      T x0  = zb0 - pmul_tan4(x3);

      zb3      = -zb3;
      zb1     -= pmul_tan4(zb3);
      T x2  = zb3 + pmul_sin4(zb1);
      T x1  = zb1 - pmul_tan4(x2);
      // Finally, the output bufferfly.
      x4       = -x4;
      x0      -= pmul_tan4(x4);
      x4      += pmul_sin4(x0);
      x0      -= pmul_tan4(x4);
      dp[0]    = x0;
      dp[7]    = x4;

      x5       = -x5;
      x1      -= pmul_tan4(x5);
      x5      += pmul_sin4(x1);
      x1      -= pmul_tan4(x5);
      dp[1]    = x1;
      dp[6]    = x5;

      x6       = -x6;
      x2      -= pmul_tan4(x6);
      x6      += pmul_sin4(x2);
      x2      -= pmul_tan4(x6);
      dp[2]    = x2;
      dp[5]    = x6;

      x7       = -x7;
      x3      -= pmul_tan4(x7);
      x7      += pmul_sin4(x3);
      x3      -= pmul_tan4(x7);
      dp[3]    = x3;
      dp[4]    = x7;
      //
    }
    //
    // Finally, loop over the columns
    dpend = target + 8;
    for(dp = target;dp < dpend;dp++) {
      // Inverse output permutation.
      T z20    =  dp[0 << 3];
      T zc0    =  dp[1 << 3];
      T z21    =  dp[2 << 3];
      T z1     = -dp[3 << 3];
      T z10    = -dp[4 << 3];
      T x45    =  dp[5 << 3];
      T z11    = -dp[6 << 3];
      T zc2    =  dp[7 << 3];
      // rotation by 45 degrees of x45,x46 to zc1,zc3.
      T z0     = x45 - pmul_tan4(z1);
      T zc3    = z1  + pmul_sin4(z0);
      T zc1    = z0  - pmul_tan4(zc3);
      // Next rotation pair.
      T z00    = z20 - pmul_tan4(z10);
      T z01    = z21 - pmul_tan2(z11);
      T zb1    = z10 + pmul_sin4(z00);
      T zb3    = z11 + pmul_sin2(z01);
      T zb0    = z00 - pmul_tan4(zb1);
      T zb2    = z01 - pmul_tan2(zb3);
      // Small butterfly. This is not exactly inverted, but inverted up to a scale factor of 2.
      zc1      = -zc1;
      zc0     -= pmul_tan4(zc1);
      z21      = zc1 + pmul_sin4(zc0);
      z20      = zc0 - pmul_tan4(z21);
      zc2      = -zc2;
      zc3     -= pmul_tan4(zc2);
      z10      = zc2 + pmul_sin4(zc3);
      z11      = zc3 - pmul_tan4(z10);
      // Rotation by 3Pi/16 and 1Pi/16.
      z00      = z20 - pmul_tan1(z10);
      z01      = z21 - pmul_tan3(z11);
      T x7     = z10 + pmul_sin1(z00);
      T x6     = z11 + pmul_sin3(z01);
      T x4     = z00 - pmul_tan1(x7);
      T x5     = z01 - pmul_tan3(x6);
      // Small butterfly again to compute x0,x1,x2,x3, again only inverted up to a scale factor of 2.
      zb2      = -zb2;
      zb0     -= pmul_tan4(zb2);
      T x3     = zb2 + pmul_sin4(zb0);
      T x0     = zb0 - pmul_tan4(x3);

      zb3      = -zb3;
      zb1     -= pmul_tan4(zb3);
      T x2     = zb3 + pmul_sin4(zb1);
      T x1     = zb1 - pmul_tan4(x2);
      // Finally, the output bufferfly.
      x4       = -x4;
      x0      -= pmul_tan4(x4);
      x4      += pmul_sin4(x0);
      x0      -= pmul_tan4(x4);
      dp[0 << 3] = x0 << preshift;
      dp[7 << 3] = x4 << preshift;

      x5       = -x5;
      x1      -= pmul_tan4(x5);
      x5      += pmul_sin4(x1);
      x1      -= pmul_tan4(x5);
      dp[1 << 3] = x1 << preshift;
      dp[6 << 3] = x5 << preshift;

      x6       = -x6;
      x2      -= pmul_tan4(x6);
      x6      += pmul_sin4(x2);
      x2      -= pmul_tan4(x6);
      dp[2 << 3] = x2 << preshift;
      dp[5 << 3] = x6 << preshift;

      x7       = -x7;
      x3      -= pmul_tan4(x7);
      x7      += pmul_sin4(x3);
      x3      -= pmul_tan4(x7);
      dp[3 << 3] = x3 << preshift;
      dp[4 << 3] = x7 << preshift;
    }
}
///
