#include "lzma_compress.hh"

#include <cstdlib>
#include <stdexcept>
#include <lzma.h>

// example from
// http://ptspts.blogspot.com/2011/11/how-to-simply-compress-c-string-with.html

std::string lzma_compress(const std::string& in, int level) {
  std::string result;
  result.resize(in.size() + (in.size() >> 2) + 128);
  size_t out_pos = 0;
  if (LZMA_OK != lzma_easy_buffer_encode(
      level, LZMA_CHECK_CRC32, NULL,
      reinterpret_cast<uint8_t*>(const_cast<char*>(in.data())), in.size(),
      reinterpret_cast<uint8_t*>(&result[0]), &out_pos, result.size()))
    throw std::runtime_error("lzma compression error");
  result.resize(out_pos);
  return result;
}
