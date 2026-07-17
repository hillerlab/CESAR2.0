/**                                                     
 * Matrix Test
 * Copyright 2017 MPI-CBG/MPI-PKS Peter Schwede
 */
#include <criterion/criterion.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "Logodd.h"

#include "Matrix.h"

struct LogoddMatrix* TestMatrix__dummy;
struct PathMatrix* TestMatrix__puppet;

Test(matrix, LogoddMatrix__str) {
  LOGODD_T dephault = .5;
  TestMatrix__dummy = LogoddMatrix__create(3, 3, dephault);

  char buffer[1024] = "";
  LogoddMatrix__str(TestMatrix__dummy, buffer);

  char cmp[1024] = "0\t+5.000000E-01\t+5.000000E-01\t+5.000000E-01\t\n1\t+5.000000E-01\t+5.000000E-01\t+5.000000E-01\t\n2\t+5.000000E-01\t+5.000000E-01\t+5.000000E-01\t\n";

  cr_assert_str_eq(buffer, cmp);
}

Test(matrix, PathMatrix__str) {
  PATH_ENTRY_T dephault = 10;
  TestMatrix__puppet = PathMatrix__create(3, 3, dephault);

  PathMatrix__set(TestMatrix__puppet, 0, 0, 3);
  PathMatrix__set(TestMatrix__puppet, 0, 1, PATH_EMPTY);
  PathMatrix__set(TestMatrix__puppet, 2, 2, 14);

  char buffer[1024] = "";
  PathMatrix__str(TestMatrix__puppet, buffer);

  char cmp[1024] = "0\t3\t10\t10\t\n1\tNA\t10\t10\t\n2\t10\t10\t14\t\n";

  cr_assert_str_eq(buffer, cmp);
  cr_assert_eq(PathMatrix__bytes(TestMatrix__puppet), 5);
}
