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


#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "denoise.c"
#include "nestim.c"

#include "lv2/lv2plug.in/ns/lv2core/lv2.h"

using namespace std;

#define NREPEL_URI "https://github.com/lucianodato/noise-repellent"
#define DEFAULT_BUFFER 4096

///---------------------------------------------------------------------

//AUXILIARY STUFF

// Works on little-endian machines only
static inline bool 
is_nan(float& value ) {
    if (((*(uint32_t *) &value) & 0x7fffffff) > 0x7f800000) {
      return true;
    }
    return false;
}

// Force already-denormal float value to zero
static inline void 
sanitize_denormal(float& value) {
    if (is_nan(value)) {
        value = 0.f;
    }
}

static inline int 
sign(float x) {
        return (x >= 0.f ? 1 : -1);
}

static inline float 
from_dB(float gdb) {
        return (exp(gdb/20.f*log(10.f)));
}

static inline float
to_dB(float g) {
        return (20.f*log10(g));
}

///---------------------------------------------------------------------

//LV2 CODE

typedef enum {
	NREPEL_INPUT  = 0,
	NREPEL_OUTPUT = 1,
	NREPEL_CAPTURE = 2,
	NREPEL_AMOUNT = 3,
	NREPEL_BUFFER = 4,

} PortIndex;

typedef struct {
	const float* input;
	float* output;

	float srate;
	
	const int* captstate;
	const float* amountreduc;
	int* bufsize;
	
} Nrepel;

static LV2_Handle
instantiate(const LV2_Descriptor*     descriptor,
            double                    rate,
            const char*               bundle_path,
            const LV2_Feature* const* features)
{
	Nrepel* nrepel = (Nrepel*)malloc(sizeof(Nrepel));
	
	nrepel->srate = rate;
			
	return (LV2_Handle)nrepel;
}

static void
connect_port(LV2_Handle instance,
             uint32_t   port,
             void*      data)
{
	Nrepel* nrepel = (Nrepel*)instance;

	switch ((PortIndex)port) {
	case NREPEL_INPUT:
		nrepel->input = (const float*)data;
		break;
	case NREPEL_OUTPUT:
		nrepel->output = (float*)data;
		break;
	case NREPEL_CAPTURE:
		nrepel->captstate = (const int*)data;
		break;
	case NREPEL_AMOUNT:
		nrepel->amountreduc = (const float*)data;
		break;
	case NREPEL_BUFFER:
		nrepel->bufsize = (int*)data;
		break;
	}
}

static void
activate(LV2_Handle instance)
{
}

static void
run(LV2_Handle instance, uint32_t n_samples)
{
	Nrepel* nrepel = (Nrepel*)instance;

	const float* input = nrepel->input;
	float* const output = nrepel->output;

	for (uint32_t pos = 0; pos < n_samples; pos++) {

			
			//finally
			output[pos] = input[pos];

	}
}

static void
deactivate(LV2_Handle instance)
{
}

static void
cleanup(LV2_Handle instance)
{
	free(instance);
}

const void*
extension_data(const char* uri)
{
	return NULL;
}

static const 
LV2_Descriptor descriptor = {
	NREPEL_URI,
	instantiate,
	connect_port,
	activate,
	run,
	deactivate,
	cleanup,
	extension_data
};

LV2_SYMBOL_EXPORT
const LV2_Descriptor*
lv2_descriptor(uint32_t index)
{
	switch (index) {
	case 0:
		return &descriptor;
	default:
		return NULL;
	}
}
