#include "fr.hpp"

using namespace barretenberg;

#define WASM_EXPORT __attribute__((visibility("default")))

extern "C" {

WASM_EXPORT void do_n_muls(size_t num, uint8_t* result)
{
    fr x = fr::one();
    for (size_t i = 0; i < num; ++i) {
        x = x * (2+i);
    }
    fr::serialize_to_buffer(x, result);
}
}