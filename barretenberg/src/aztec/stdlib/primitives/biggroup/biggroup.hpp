#pragma once

#include "../bigfield/bigfield.hpp"
#include "../field/field.hpp"

namespace plonk {
namespace stdlib {

// ( ͡° ͜ʖ ͡°)
template <typename Composer, class Fq, class Fr, class Params> class element {
  public:
    element();
    element(const Fq& x, const Fq& y);

    element(const element& other);
    element(element&& other);

    element& operator=(const element& other);
    element& operator=(element&& other);

    element operator+(const element& other) const;
    element operator-(const element& other) const;
    element operator*(const Fr& other) const;

    element dbl() const;
    element montgomery_ladder(const element& other) const;

    static element twin_mul(const element& base_a, const Fr& scalar_a, const element& base_b, const Fr& scalar_b);

    static element quad_mul(const element& base_a,
                            const Fr& scalar_a,
                            const element& base_b,
                            const Fr& scalar_b,
                            const element& base_c,
                            const Fr& scalar_c,
                            const element& base_d,
                            const Fr& scalar_d);

    static std::vector<bool_t<Composer>> compute_naf(const Fr& scalar);
    Composer* get_context(const element& other) const
    {
        if (x.context != nullptr) {
            return x.context;
        }
        if (y.context != nullptr) {
            return y.context;
        }
        if (other.x.context != nullptr) {
            return other.x.context;
        }
        if (other.y.context != nullptr) {
            return other.y.context;
        }
    }

    Fq x;
    Fq y;
};
} // namespace stdlib
} // namespace plonk

#include "biggroup_impl.hpp"