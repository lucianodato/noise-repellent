#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "lv2/lv2plug.in/ns/lv2core/lv2.h"

#define NREPEL_URI "https://github.com/lucianodato/noise-repellent"

typedef enum {
	NREPEL_INPUT  = 0,
	NREPEL_OUTPUT = 1,

} PortIndex;

typedef struct {
	float* input;
	float* output;

	float srate;

} Nrepel;



static LV2_Handle
instantiate(const LV2_Descriptor*     descriptor,
            double                    rate,
            const char*               bundle_path,
            const LV2_Feature* const* features)
{
	int i;
	Nrepel* nrepel = (Nreduc*)malloc(sizeof(Nrepel));
	
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
		zameq2->input = (float*)data;
		break;
	case NREPEL_OUTPUT:
		zameq2->output = (float*)data;
		break;
	}
}

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

static void
activate(LV2_Handle instance)
{
}

static void
run(LV2_Handle instance, uint32_t n_samples)
{
	Nrepel* nrepel = (Nrepel*)instance;

	const float* const input  = nrepel->input;
	float* const       output = nrepel->output;

	for (uint32_t pos = 0; pos < n_samples; pos++) {
	
			float in = input[pos];

			sanitize_denormal(in);

			//finally
			output[pos]=in;

		}
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

static const LV2_Descriptor descriptor = {
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
