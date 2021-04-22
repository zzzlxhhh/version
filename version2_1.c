int DSPF_sp_qrd(const int Nrows, const int Ncols, float *A, float *Q, float *R,
                float *u) {

  int row, col, i, k, loop_count, num, vBlock_index, trans_index;
  float alpha, scale, scale_tmp, sum, sump, norm_sqr, diag;
  num = Nrows * Ncols;

  int Nvecs = Nrows / 16; // R的向量数
  int Ri = -1;            // Vec matrix_index
  vector float tmp1_1,tmp1_2,tmp1_3;
  vector float tmp2_1,tmp2_2,tmp2_3;
  vector float tmp3_1,tmp3_2,tmp3_3;
  vector float tmp4_1,tmp4_2,tmp4_3;
  vector float tmp5_1,tmp5_2,tmp5_3;
  vector float tmp6_1,tmp6_2,tmp6_3;

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
        /* R=Q1*R */
        for (i = col + 1; i < Ncols; i++) {
          /*sum = com_norm_sqr(&r[i * Nvecs + Ri], &uv[Ri], Nvecs - Ri,
                             vBlock_index);*/
        

          sum *= scale;
          v_mul_sub(&r[i * Nvecs + Ri], &uv[Ri], sum, Nvecs - Ri, vBlock_index);
        }
        /* Q=A*Q1 */
        for (i = 0; i < Nrows; i++) {
          sum = com_norm_sqr(&q[i * Nvecs + Ri], &uv[Ri], Nvecs - Ri,
                             vBlock_index);
          sum *= scale;
          v_mul_sub(&q[i * Nvecs + Ri], &uv[Ri], sum, Nvecs - Ri, vBlock_index);
        }
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
