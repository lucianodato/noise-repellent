/*
noise-repellent -- Noise Reduction LV2

Copyright 2016 Luciano Dato <lucianodato@gmail.com>

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

/**
* \file circular_buffer.c
* \author Luciano Dato
* \brief Contains a very basic circular buffer abstraction
*/

#include <malloc.h>
#include <float.h>

typedef struct
{
  int write_pointer;
  int read_pointer;
  int size;
  float* data;
} circular_buffer;

static circular_buffer*
cb_init(int size)
{
  circular_buffer *cb = (circular_buffer*)malloc(sizeof(circular_buffer));
  cb->size = size;
  cb->write_pointer = 0;
  cb->read_pointer = 0;
  cb->data = (float*)calloc(cb->size, sizeof(float));

  return cb;
}

static void
cb_free(circular_buffer *cb)
{
  free(cb->data);
  free(cb);
}

static void
cb_reset(circular_buffer *cb)
{
  cb->read_pointer = cb->write_pointer;
}

static bool
cb_is_full(circular_buffer *cb)
{
  return((cb->write_pointer + 1) % cb->size == cb->read_pointer );
}

static bool
cb_is_empty(circular_buffer *cb)
{
  return(cb->read_pointer == cb->write_pointer);
}

static bool
cb_write_one(circular_buffer *cb, float value)
{
  bool is_full = cb_is_full(cb);
  cb->data[cb->write_pointer] = value;
  cb->write_pointer++;
  cb->write_pointer %= cb->size;
  return is_full;
}

static bool
cb_read_one(circular_buffer *cb, float *value)
{
  bool is_empty = cb_is_empty(cb);
  *value = cb->data[cb->read_pointer];
  cb->read_pointer++;
  cb->read_pointer %= cb->size;
  return(is_empty);
}
