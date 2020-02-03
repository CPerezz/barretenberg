#include "./lagrange_base.hpp"

namespace barretenberg
{
namespace lagrange_base
{
    void transform_srs(g1::affine_element* monomials, g1::affine_element* lagrange_base, const size_t degree)
    {
        fr::field_t n_inv = fr::to_montgomery_form({{ static_cast<uint64_t>(degree), 0, 0, 0 }});
        n_inv = fr::invert(n_inv);
        for (size_t i = 0; i < degree; ++i)
        {
            lagrange_base[i] = g1::group_exponentiation(monomials[0], n_inv);
        }
    }
}
}