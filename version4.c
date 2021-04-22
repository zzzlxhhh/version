#include <math.h>
#include <string.h>

#include "DSPF_sp_qrd_cn.h"

/*
循环展开,计算访存解耦
*/
void part_op(int vBlock_index) {
  switch (vBlock_index) {
  case 0: {
    mov_to_vlr(0xFFFF);
    break;
  }
  case 1: {
    mov_to_vlr(0xFFFE);
    break;
  }
  case 2: {
    mov_to_vlr(0xFFFC);
    break;
  }
  case 3: {
    mov_to_vlr(0xFFF8);
    break;
  }
  case 4: {
    mov_to_vlr(0xFFF0);
    break;
  }
  case 5: {
    mov_to_vlr(0xFFE0);
    break;
  }
  case 6: {
    mov_to_vlr(0xFFC0);
    break;
  }
  case 7: {
    mov_to_vlr(0xFF80);
    break;
  }
  case 8: {
    mov_to_vlr(0xFF00);
    break;
  }
  case 9: {
    mov_to_vlr(0xFE00);
    break;
  }
  case 10: {
    mov_to_vlr(0xFC00);
    break;
  }
  case 11: {
    mov_to_vlr(0xF800);
    break;
  }
  case 12: {
    mov_to_vlr(0xF000);
    break;
  }
  case 13: {
    mov_to_vlr(0xE000);
    break;
  }
  case 14: {
    mov_to_vlr(0xC000);
    break;
  }
  case 15: {
    mov_to_vlr(0x8000);
    break;
  }
  }
}
void part_op_uv(int vBlock_index) {
  switch (vBlock_index) {
  case 15: {
    mov_to_vlr(0x8000);
    break;
  }
  case 14: {
    mov_to_vlr(0x4000);
    break;
  }
  case 13: {
    mov_to_vlr(0x2000);
    break;
  }
  case 12: {
    mov_to_vlr(0x1000);
    break;
  }
  case 11: {
    mov_to_vlr(0x0800);
    break;
  }
  case 10: {
    mov_to_vlr(0x0400);
    break;
  }
  case 9: {
    mov_to_vlr(0x0200);
    break;
  }
  case 8: {
    mov_to_vlr(0x0100);
    break;
  }
  case 7: {
    mov_to_vlr(0x0080);
    break;
  }
  case 6: {
    mov_to_vlr(0x0040);
    break;
  }
  case 5: {
    mov_to_vlr(0x0020);
    break;
  }
  case 4: {
    mov_to_vlr(0x0010);
    break;
  }
  case 3: {
    mov_to_vlr(0x0008);
    break;
  }
  case 2: {
    mov_to_vlr(0x0004);
    break;
  }
  case 1: {
    mov_to_vlr(0x0002);
    break;
  }
  case 0: {
    mov_to_vlr(0x0001);
    break;
  }
  }
}
//更新uv、R主对角元?  也可尝试用svr更新
void update_uv_R(vector float *v, vector float *uv, int vBlock_index,
                 float alpha) {
  vector float vtmp;
  vtmp = vec_svbcast(alpha);
  part_op_uv(vBlock_index);

  // update uv u[col] = R[col + col * Ncols] + alpha
  uv[0] = vec_add(uv[0], vtmp);

  // update R   R[col + col * Ncols] = -alpha;
  v[0] = vec_sub(vtmp, v[0]);

  mov_to_vlr(0xFFFF);
}

float reduce_16(vector float *v) {
  float array[16];
  M7002_datatrans(v, array, 64);
  float res = 0;
  res = array[0] + array[1] + array[2] + array[3] + array[4] + array[5] +
        array[6] + array[7] + array[8] + array[9] + array[10] + array[11] +
        array[12] + array[13] + array[14] + array[15];

  return res;
}

//计算第col列的模长，顺便更新一部分uv和R的�?
float norm2(vector float *v, vector float *uv, int vBlocks, int vBlock_index) {
  int i, loop, b1, b2;
  float res;

  vector float vtmp, vzero, vtmp1, vtmp2;
  vtmp = vec_svbcast(0.0f);
  vzero = vtmp;
  vtmp1 = vtmp;
  vtmp2 = vtmp;

  part_op(vBlock_index);
  uv[0] = vec_muli(uv[0], vzero);
  uv[0] = vec_add(uv[0], v[0]);

  vtmp = vec_mula(v[0], v[0], vtmp);
  v[0] = vec_muli(v[0], vzero);
  mov_to_vlr(0xFFFF);

  loop = (vBlocks - 1) % 3;
  loop = vBlocks - loop;
  i = 1;
  for (; i < loop; i++) {
    b1 = ++i;
    b2 = ++i;

    vtmp = vec_mula(v[i], v[i], vtmp);
    vtmp1 = vec_mula(v[b1], v[b1], vtmp1);
    vtmp2 = vec_mula(v[b2], v[b2], vtmp2);

    uv[i] = v[i];
    uv[b1] = v[b1];
    uv[b2] = v[b2];

    v[i] = vzero;
    v[b1] = vzero;
    v[b2] = vzero;
  }
  for (; i < vBlocks; i++) {
    vtmp = vec_mula(v[i], v[i], vtmp);
    uv[i] = v[i];
    v[i] = vzero;
  }
  vtmp = vec_add(vtmp, vtmp1);
  vtmp = vec_add(vtmp, vtmp2);
  res = reduce_16(&vtmp);
  return res;
}
void update_qr(vector float *v1, vector float *v2, int vBlocks,
               int vBlock_index, int rows, int Nvecs, float scale) {
  int i, j, loop, in_loop, index;
  vector float sum1, sum2, sum3, vzero;
  vzero = vec_svbcast(0.0f);
  float res1, res2, res3;

  loop = rows - rows % 3;
  in_loop = (vBlocks - 1) % 3;
  in_loop = vBlocks - in_loop;

  for (i = 0; i < loop; i += 3) {
    sum1 = vzero;
    sum2 = sum1;
    sum3 = sum1;

    index = i * Nvecs;
    part_op(vBlock_index);
    sum1 = vec_mula(v1[index], v2[0], sum1);
    sum2 = vec_mula(v1[index + Nvecs], v2[0], sum2);
    sum3 = vec_mula(v1[index + 2 * Nvecs], v2[0], sum3);
    mov_to_vlr(0xFFFF);
    for (j = 1; j < in_loop; j+=3) {
      index = j + i * Nvecs;
      sum1 = vec_mula(v1[index], v2[j], sum1);
      sum2 = vec_mula(v1[index + Nvecs], v2[j], sum2);
      sum3 = vec_mula(v1[index + 2 * Nvecs], v2[j], sum3);
      index++;
      sum1 = vec_mula(v1[index], v2[j+1], sum1);
      sum2 = vec_mula(v1[index + Nvecs], v2[j+1], sum2);
      sum3 = vec_mula(v1[index + 2 * Nvecs], v2[j+1], sum3);
      index++;
      sum1 = vec_mula(v1[index], v2[j+2], sum1);
      sum2 = vec_mula(v1[index + Nvecs], v2[j+2], sum2);
      sum3 = vec_mula(v1[index + 2 * Nvecs], v2[j+2], sum3);
    }
    for (; j < vBlocks; j++) {
      index = j + i * Nvecs;
      sum1 = vec_mula(v1[index], v2[j], sum1);
      sum2 = vec_mula(v1[index + Nvecs], v2[j], sum2);
      sum3 = vec_mula(v1[index + 2 * Nvecs], v2[j], sum3);
    }
    res1 = reduce_16(&sum1);
    res2 = reduce_16(&sum2);
    res3 = reduce_16(&sum3);

    res1 *= scale;
    sum1 = vec_svbcast(res1);
    res2 *= scale;
    sum2 = vec_svbcast(res2);
    res3 *= scale;
    sum3 = vec_svbcast(res3);

    index = i * Nvecs;
    part_op(vBlock_index);
    v1[index] = vec_mulb(v2[0], sum1, v1[index]);
    v1[index] = vec_neg(v1[index]);
    v1[index + Nvecs] = vec_mulb(v2[0], sum2, v1[index + Nvecs]);
    v1[index + Nvecs] = vec_neg(v1[index + Nvecs]);
    v1[index + Nvecs * 2] =vec_mulb(v2[0], sum3, v1[index + Nvecs * 2]);
    v1[index + Nvecs * 2] = vec_neg(v1[index + Nvecs * 2]);
    mov_to_vlr(0xFFFF);
    for (j = 1; j < in_loop; j+=3) {
      index = j + i * Nvecs;
      v1[index] = vec_mulb(v2[j], sum1, v1[index]);
      v1[index] = vec_neg(v1[index]);
      v1[index + Nvecs] = vec_mulb(v2[j], sum2, v1[index + Nvecs]);
      v1[index + Nvecs] = vec_neg(v1[index + Nvecs]);
      v1[index + Nvecs * 2] =vec_mulb(v2[j], sum3, v1[index + Nvecs * 2]);
      v1[index + Nvecs * 2] = vec_neg(v1[index + Nvecs * 2]);
      index++;
      v1[index] = vec_mulb(v2[j+1], sum1, v1[index]);
      v1[index] = vec_neg(v1[index]);
      v1[index + Nvecs] = vec_mulb(v2[j+1], sum2, v1[index + Nvecs]);
      v1[index + Nvecs] = vec_neg(v1[index + Nvecs]);
      v1[index + Nvecs * 2] =vec_mulb(v2[j+1], sum3, v1[index + Nvecs * 2]);
      v1[index + Nvecs * 2] = vec_neg(v1[index + Nvecs * 2]);
      index++;
      v1[index] = vec_mulb(v2[j+2], sum1, v1[index]);
      v1[index] = vec_neg(v1[index]);
      v1[index + Nvecs] = vec_mulb(v2[j+2], sum2, v1[index + Nvecs]);
      v1[index + Nvecs] = vec_neg(v1[index + Nvecs]);
      v1[index + Nvecs * 2] =vec_mulb(v2[j+2], sum3, v1[index + Nvecs * 2]);
      v1[index + Nvecs * 2] = vec_neg(v1[index + Nvecs * 2]);
    }
    for (; j < vBlocks; j++) {
      index = j + i * Nvecs;
      v1[index] = vec_mulb(v2[j], sum1, v1[index]);
      v1[index] = vec_neg(v1[index]);
      v1[index + Nvecs] = vec_mulb(v2[j], sum2, v1[index + Nvecs]);
      v1[index + Nvecs] = vec_neg(v1[index + Nvecs]);
      v1[index + Nvecs * 2] =vec_mulb(v2[j], sum3, v1[index + Nvecs * 2]);
      v1[index + Nvecs * 2] = vec_neg(v1[index + Nvecs * 2]);
    }
  }
  //余下不能三行操作部分
  for (; i < rows; i++) {
    sum1 = vzero;

    index = i * Nvecs;
    part_op(vBlock_index);
    sum1 = vec_mula(v1[index], v2[0], sum1);
    mov_to_vlr(0xFFFF);

    for (j = 1; j < in_loop; j += 3) {
      index = j + i * Nvecs;
      sum1 = vec_mula(v1[index], v2[j], sum1);
      sum1 = vec_mula(v1[index + 1], v2[j], sum1);
      sum1 = vec_mula(v1[index + 2], v2[j], sum1);
    }
    for (; j < vBlocks; j++) {
      index = j + i * Nvecs;
      sum1 = vec_mula(v1[index], v2[j], sum1);
    }
    res1 = reduce_16(&sum1);
    res1 *= scale;
    sum1 = vec_svbcast(res1);

    index = i * Nvecs;
    part_op(vBlock_index);
    v1[index] = vec_mulb(v2[0], sum1, v1[index]);
    v1[index] = vec_neg(v1[index]);
    mov_to_vlr(0xFFFF);
    for (j = 1; j < in_loop; j+=3) {
      index = j + i * Nvecs;
      v1[index] = vec_mulb(v2[j], sum1, v1[index]);
      v1[index] = vec_neg(v1[index]);
      index++;
      v1[index] = vec_mulb(v2[j+1], sum1, v1[index]);
      v1[index] = vec_neg(v1[index]);
      index++;
      v1[index] = vec_mulb(v2[j+2], sum1, v1[index]);
      v1[index] = vec_neg(v1[index]);
    }
    for (; j < vBlocks; j++) {
      index = j + i * Nvecs;
      v1[index] = vec_mulb(v2[j], sum1, v1[index]);
      v1[index] = vec_neg(v1[index]);
    }
  }
}

//处理16对齐
int DSPF_sp_qrd(const int Nrows, const int Ncols, float *A, float *Q, float *R,
                float *u) {

  int row, col, i, k, loop_count, num, vBlock_index, trans_index;
  float alpha, scale, scale_tmp, sum, sump, norm_sqr, diag;
  num = Nrows * Ncols;

  int Nvecs = Nrows / 16; // R的向量数
  int Ri = -1;            // Vec matrix_index

  memset(Q, 0.0, sizeof(float) * Nrows * Nrows);
  for (row = 0; row < Nrows; row++) {
    Q[row + row * Nrows] = 1.0;
  }

  vector float *r = (vector float *)vmalloc(num * 4);
  vector float *uv = (vector float *)vmalloc(Nrows * 4);
  vector float *q = (vector float *)vmalloc(Nrows * Nrows * 4);
  vector float castv;
  //转置传输
  trans_index = Ncols * 4 - 4;
  for (col = 0; col < Ncols; col++) {
    M7002_datatrans_index(&A[col], &r[col * Nvecs], Nrows, 1, trans_index);
  }
  M7002_datatrans(Q, q, Nrows * Nrows * 4);

  if (Nrows <= Ncols) {
    loop_count = Nrows - 2;
  } else {
    loop_count = Ncols - 1;
  }
  for (col = 0; col <= loop_count; col++) {
    /*
    if (col % 16 == 0)
        Ri++; // 在向量矩阵中的行index
        */
    Ri = (int)col * 0.0625;
    sum = 0;
    vBlock_index = col % 16;
    sum = norm2(&r[col * Nvecs + Ri], &uv[Ri], Nvecs - Ri, vBlock_index);
    if (sum != 0) {
      vmemcpy(&diag, (float *)(&uv[Ri]) + vBlock_index, 4);
      // M7002_datatrans((float*)(&uv[Ri])+vBlock_index,&diag,4);
      // diag=darray[vBlock_index];
      alpha = sqrt(sum);
      if (diag >= 0) {
        alpha = -alpha;
      }

      update_uv_R(&r[col * Nvecs + Ri], &uv[Ri], vBlock_index, alpha);
      // norm_sqr = com_norm_sqr(uv + Ri, uv + Ri, Nvecs - Ri, vBlock_index);

      // if (alpha * u[col] != 0.0)
      scale_tmp = alpha * (diag + alpha);

      if (scale_tmp != 0.0) {
        scale = 1 / scale_tmp;
        update_qr(&r[(col + 1) * Nvecs + Ri], &uv[Ri], Nvecs - Ri, vBlock_index,
                  Ncols - col - 1, Nvecs, scale);
        update_qr(&q[Ri], &uv[Ri], Nvecs - Ri, vBlock_index, Nrows, Nvecs,
                  scale);
      } /* if (norm_sqr!=0) */
    }   /* if (sum!=0) */

  } /* for (col=0;col<=loop_count;col++) */
  // M7002_mat_transpose(r, R, Ncols,Nrows, 0);
  trans_index = trans_index << 16;
  for (col = 0; col < Ncols; col++) {
    M7002_datatrans_index(&r[col * Nvecs], &R[col], Nrows, 1, trans_index);
  }
  M7002_datatrans(q, Q, Nrows * Nrows * 4);
  M7002_datatrans(uv, u, Nrows * 4);
  vfree(r);
  vfree(uv);
  vfree(q);
  return 0;
}
