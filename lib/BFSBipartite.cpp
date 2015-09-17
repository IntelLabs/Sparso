#include "BFSBipartite.h"

using namespace SpMP;

extern "C" {

void bfsBipartite(CSR& A, CSR& AT, int *rowPerm, int *rowInversePerm, int *colPerm, int *colInversePerm)
{
  int base = A.getBase();

  bool *visited1 = new bool[A.m];
  bool *visited2 = new bool[A.n];

  for (int i = 0; i < A.m; ++i) {
    visited1[i] = false;
  }
  for (int i = 0; i < A.n; ++i) {
    visited2[i] = false;
  }

  int rowPermIdx = 0;
  int colPermIdx = 0;

  for (int i = 0; i < A.m; ++i) {
    if (!visited1[i]) {
      std::vector<int> q1, q2;
      q1.push_back(i);
      visited1[i] = true;

      while (!q1.empty()) {
        for (int u : q1) {
          rowPerm[u] = rowPermIdx;
          rowInversePerm[rowPermIdx++] = u;
          for (int j = A.rowptr[u] - base; j < A.rowptr[u + 1] - base; ++j) {
            int v = A.colidx[j] - base;
            if (!visited2[v]) {
              visited2[v] = true;
              q2.push_back(v);
            }
          }
        }
        q1.clear();

        if (q2.empty()) {
          break;
        }

        for (int u : q2) {
          colPerm[u] = colPermIdx;
          colInversePerm[colPermIdx++] = u;
          for (int j = AT.rowptr[u] - base; j < AT.rowptr[u + 1] - base; ++j) {
            int v = AT.colidx[j] - base;
            if (!visited1[v]) {
              visited1[v] = true;
              q1.push_back(v);
            }
          }
        }
        q2.clear();
      } // while q is not empty
    } // for each connected component
  }
  for (int i = 0; i < A.n; ++i) {
    if (!visited2[i]) {
      std::vector<int> q1, q2;
      q2.push_back(i);
      visited2[i] = true;

      while (!q2.empty()) {
        for (int u : q2) {
          colPerm[u] = colPermIdx;
          colInversePerm[colPermIdx++] = u;
          for (int j = AT.rowptr[u] - base; j < AT.rowptr[u + 1] - base; ++j) {
            int v = AT.colidx[j] - base;
            if (!visited1[v]) {
              visited1[v] = true;
              q1.push_back(v);
            }
          }
        }
        q2.clear();

        if (q1.empty()) {
          break;
        }

        for (int u : q1) {
          rowPerm[u] = rowPermIdx;
          colInversePerm[rowPermIdx++] = u;
          for (int j = A.rowptr[u] - base; j < A.rowptr[u + 1] - base; ++j) {
            int v = A.colidx[j] - base;
            if (!visited2[v]) {
              visited2[v] = true;
              q2.push_back(v);
            }
          }
        }
        q1.clear();
      } // while q is not empty
    } // for each connected component
  }

  assert(rowPermIdx == A.m);
  assert(colPermIdx == A.n);
}

} // extern "C"
