#pragma once
#include <algorithm>
#include <cmath>
#include <vector>

#ifndef THE_MATRIX_DEFINED
using the_matrix_type = int;
using the_matrix = std::vector<std::vector<the_matrix_type>>;
#define THE_MATRIX_DEFINED
#endif

template <typename T>
void DCT_bin(T &a0, T &a1, T &a2, T &a3, T &a4, T &a5, T &a6, T &a7) {
  T tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp10, tmp11, tmp12, tmp13;
  tmp7 = a0 - a7;
  tmp0 = a0 - ((tmp7) / 2);

  tmp6 = a1 - a6;
  tmp1 = a1 - ((tmp6) / 2);

  tmp5 = a2 - a5;
  tmp2 = a2 - ((tmp5) / 2);

  tmp4 = a3 - a4;
  tmp3 = a3 - ((tmp4) / 2);

  /* Even part */
  tmp13 = tmp0 - tmp3;
  tmp10 = tmp0 - ((tmp13) / 2);

  tmp12 = tmp1 - tmp2;
  tmp11 = tmp1 - ((tmp12) / 2);

  a4 = tmp10 - tmp11;
  a0 = tmp10 - ((a4) / 2);

  /*1/2, -1/2: alter the sign to get positive scaling factor */
  /* new version: 7/16 and -3/8 */
  a6 = ((tmp13) / 2) - tmp12;
  a2 = tmp13 - ((a6) / 2);

  /* Odd part */

  /* pi/4 = -1/2u 3/4d -1/2u*/
  tmp10 = tmp5 - ((tmp6) / 2);
  tmp6 = tmp6 + tmp10 - ((tmp10) / 4);
  tmp5 = ((tmp6) / 2) - tmp10;

  // butterflies:
  tmp10 = tmp4 + tmp5;
  tmp11 = ((tmp10) / 2) - tmp5;

  tmp13 = tmp6 + tmp7;
  tmp12 = ((tmp13) / 2) - tmp6;

  /* 7pi/16 = 1/4u -1/4d: alter the sign to get positive scaling factor */
  a7 = ((tmp13) / 4) - tmp10;
  a1 = tmp13 - ((a7) / 4);

  /* 3pi/16 = */
  /* new version: 1, -1/2 */
  a5 = tmp11 + tmp12;
  a3 = tmp12 - ((a5) / 2);
}

template <typename T>
void inv_DCT_bin(T &a0, T &a1, T &a2, T &a3, T &a4, T &a5, T &a6, T &a7) {
  T tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp10, tmp11, tmp12, tmp13;

  /* X[0] and X[4] */
  tmp10 = a0 + ((a4) / 2);
  tmp11 = tmp10 - a4;

  /* X[6] and X[2]:1/2, 1/2 */
  tmp13 = a2 + ((a6) / 2);
  tmp12 = ((tmp13) / 2) - a6;

  // lossless binDCT: use new nutterflies.
  tmp0 = tmp10 + ((tmp13) / 2);
  tmp3 = tmp0 - tmp13;

  tmp1 = tmp11 + ((tmp12) / 2);
  tmp2 = tmp1 - tmp12;

  /* 7pi/16 = -1/4d 1/4u */
  tmp13 = a1 + ((a7) / 4);
  tmp10 = ((tmp13) / 4) - a7;

  /* 3pi/16 = 1/2d -1u */
  tmp12 = a3 + ((a5) / 2);
  tmp11 = a5 - tmp12;

  // lossless binDCT: use new butterflies.
  tmp5 = ((tmp10) / 2) - tmp11;
  tmp4 = tmp10 - tmp5;

  tmp6 = ((tmp13) / 2) - tmp12;
  tmp7 = tmp13 - tmp6;

  /* pi/4 = -1/2u -3/4d 1/2u */
  tmp5 = ((tmp6) / 2) - tmp5;
  tmp6 = tmp6 - tmp5 + ((tmp5) / 4);
  tmp5 = tmp5 + ((tmp6) / 2);

  /* last stage: butterfly */

  /* Final output stage: scale down by a factor of 8 and range-limit */
  tmp10 = tmp0 + ((tmp7) / 2);
  tmp11 = tmp10 - tmp7;
  a0 = tmp10;
  a7 = tmp11;

  tmp10 = tmp1 + ((tmp6) / 2);
  tmp11 = (tmp10 - tmp6);
  a1 = tmp10;
  a6 = tmp11;

  tmp10 = tmp2 + ((tmp5) / 2);
  tmp11 = (tmp10 - tmp5);
  a2 = tmp10;
  a5 = tmp11;

  tmp10 = tmp3 + ((tmp4) / 2);
  tmp11 = (tmp10 - tmp4);
  a3 = tmp10;
  a4 = tmp11;
}

void DCT_bin(the_matrix &a) {

  auto height = a.size() + 7 & ~7;
  if (height == 0)
    return;
  auto width = a[0].size() + 7 & ~7;
  if (width == 0)
    return;

  for (int y = (a.resize(height, a.back()), 0); y < height; y += 1)
    for (int x = (a[y].resize(width, a[y].back()), 0); x < width; x += 8) {
      DCT_bin(a[y][x], a[y][x + 1], a[y][x + 2], a[y][x + 3], a[y][x + 4],
              a[y][x + 5], a[y][x + 6], a[y][x + 7]);
    }

  for (int y = 0; y < height; y += 8)
    for (int x = 0; x < width; x += 1) {
      DCT_bin(a[y][x], a[y + 1][x], a[y + 2][x], a[y + 3][x], a[y + 4][x],
              a[y + 5][x], a[y + 6][x], a[y + 7][x]);
    }
}
static void inv_DCT_bin(the_matrix &a) {
  const auto height = a.size();
  const auto width = a[0].size();
  for (int y = 0; y < height; y += 8)
    for (int x = 0; x < width; x += 1) {
      inv_DCT_bin(a[y][x], a[y + 1][x], a[y + 2][x], a[y + 3][x], a[y + 4][x],
                  a[y + 5][x], a[y + 6][x], a[y + 7][x]);
    }
  for (int y = 0; y < height; y += 1)
    for (int x = 0; x < width; x += 8) {
      inv_DCT_bin(a[y][x], a[y][x + 1], a[y][x + 2], a[y][x + 3], a[y][x + 4],
                  a[y][x + 5], a[y][x + 6], a[y][x + 7]);
    }
}
