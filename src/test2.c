#include <stdio.h>
#include <stdlib.h>

#define TEST(var,type) \
  printf("\033[36m" #var "\033[0m = " type "\n", var);

int main(int argc, char* argv[]) {
  size_t buf_len = 1 << 7;
  char* m = malloc(1<<10);
  size_t m_len = 1<<10;
  size_t m_used = 0;
  // FILE* pipe = popen("cat test.txt","r");
  FILE* pipe = popen("unxz -c test.dat.xz","r");
  if (!pipe) return 1;

  for (;;) {
    size_t n = fread(m+m_used, 1, buf_len, pipe);
    m_used += n;
    TEST(m_len,"%zu")
    TEST(m_used,"%zu")
    if (n < buf_len) break;
    if (m_used >= m_len)
      m = realloc(m,(m_len<<=1)+1);
  }

  pclose(pipe);
  free(m);
}

