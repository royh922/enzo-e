#include <stdio.h>
#include "grackle.h"

grackle_version get_grackle_version(void) {
  grackle_version out;
  out.version = "3.3.1-dev";
  out.branch = "main";
  out.revision = "4798687a1c0e46489af679a555e22ad3f992fede";
  return out;
}

void auto_show_flags(FILE *fp) {
  fprintf (fp, "%s\n", "\n  CC = /usr/bin/clang\n  FC = /opt/homebrew/bin/gfortran\n  LD = /usr/bin/ld\n");
}

void auto_show_config(FILE *fp) {
  fprintf (fp, "%s\n", "\n   Built with CMake\n\n   GRACKLE_USE_DOUBLE                  : ON\n   GRACKLE_USE_OPENMP                  : OFF\n   BUILD_TYPE                          : Release\n");
}
