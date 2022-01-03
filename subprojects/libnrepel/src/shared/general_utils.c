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

inline float sanitize_denormal(float value) {
  if (!isnormal(value)) {
    value = 0.F;
  }
  return value;
}

inline float from_db_to_coefficient(const float gain_db) {
  return expf(gain_db / 10.F * logf(10.F));
}