/*
noise-repellent -- Noise Reduction LV2

Copyright 2021 Luciano Dato <lucianodato@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/
*/

#include "general_utils.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>

inline float sanitize_denormal(float value) {
  if (!isnormal(value)) {
    value = 0.F;
  }
  return value;
}

inline float from_db_to_coefficient(const float value_db) {
  return expf(value_db / 10.F * logf(10.F));
}

int get_next_divisible_two(int number) {
  int q = number / 2;
  int n1 = 2 * q;
  int n2 = (number * 2) > 0 ? (2 * (q + 1)) : (2 * (q - 1));
  if (abs(number - n1) < abs(number - n2)) {
    return n1;
  }

  return n2;
}

int get_next_power_two(int number) {
  return (int)roundf(powf(2.F, ceilf(log2f((float)number))));
}