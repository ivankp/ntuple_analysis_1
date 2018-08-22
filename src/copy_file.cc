#include "copy_file.hh"

#include <sys/sendfile.h> // sendfile
#include <fcntl.h>        // open
#include <unistd.h>       // close
#include <sys/stat.h>     // fstat
#include <sys/types.h>    // fstat

void copy_file(const char* from, const char* to) {
  const int _from = open(from, O_RDONLY, 0);
  const int _to   = open(to, O_WRONLY | O_CREAT, 0644);

  struct stat stat_source;
  fstat(_from, &stat_source);

  sendfile(_to, _from, 0, stat_source.st_size);

  close(_from);
  close(_to);
}
