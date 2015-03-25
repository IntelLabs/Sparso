#include <cassert>
#include <climits>
#include <vector>
#include <algorithm>

#include "../CSR.hpp"

using namespace std;

void getInversePerm(int *inversePerm, const int *perm, int n);

void CSR::getRCMPermutationNew(int *perm, int *inversePerm, int source /*=-1*/) const
{
  // 1. Start vertex
  // TODO: pseudo diameter heuristic
  assert(source >= 0 && source < m);

  // 2. BFS
  int *levels = new int[m];
  for (int i = 0; i < m; ++i) {
    levels[i] = INT_MAX;
  }

  int level = 0;
  int nextId = 0;

  levels[source] = level;
  inversePerm[nextId++] = source;

  vector<int> q, nextQ;
  q.push_back(source);

  while (!q.empty()) {
    ++level;
    for (auto u = q.begin(); u != q.end(); ++u) {
      int index = nextQ.size();
      assert(levels[u] == level);

      for (int j = rowPtr[*u]; j < rowPtr[*u + 1]; ++j) {
        int v = colIdx[j];
        if (levels[v] > levels[*u] + 1) {
          levels[v] = levels[*u] + 1;
          nextQ.push_back(v);
        }
      }
      sort(&nextQ[index], &nextQ[nextQ.size()]);
    } // for each u in current level

    for (auto v = nextQ.begin(); v != nextQ.end(); ++v) {
      inversePerm[nextId++] = *v;
    }
    q.clear();
    q.swap(nextQ);
  } // while !q.empty()

  // 3. Reorder
  delete[] levels;

  int *temp = new int[m];
  reverse_copy(inversePerm, inversePerm + m, temp);
  copy(temp, temp + m, inversePerm);
  delete[] temp;

  getInversePerm(perm, inversePerm, n);
}
