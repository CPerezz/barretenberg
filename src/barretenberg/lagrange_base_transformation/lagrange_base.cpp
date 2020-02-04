#include "./lagrange_base.hpp"
#include <iostream>
namespace barretenberg {
namespace lagrange_base {

std::vector<g1::element> g1fft(const g1::element* monomials, size_t size, barretenberg::fr::field_t root, size_t offset)
{

        std::cout << "size:" << size << std::endl;
    std::vector<g1::element> result(size);
    if (size == 1) {
        result[0] = monomials[0];
        return result;
    }

    auto new_root = fr::mul(root, root);
    auto odd = g1fft(monomials + offset, size / 2, new_root, offset * 2);
    auto even = g1fft(monomials, size / 2, new_root, offset * 2);
        std::cout << "size:" << size << std::endl;
    auto current = root;

    for (size_t i = 0; i < size / 2; ++i) {
        fr::print(fr::from_montgomery_form(current));
        auto temp = g1::group_exponentiation(odd[i], current);
        g1::element temp2;
        g1::add(even[i], temp, temp2);
        result[i] = temp2;
        g1::__neg(temp, temp);
        g1::add(even[i], temp, temp2);
        result[size / 2 + i] = temp2;
        current = fr::mul(root, current);
    }
    return result;
}

void transform_srs(g1::affine_element* monomials, g1::affine_element* lagrange_base_affine, const size_t degree)
{
    barretenberg::evaluation_domain domain(degree);
    barretenberg::fr::field_t root = domain.root_inverse;
    std::vector<g1::element> monomials_jac;
    for (size_t i = 0; i < degree; ++i) {
        monomials_jac.push_back({ monomials[i].x, monomials[i].y, fq::one });
    }

    auto lagrange_jac = g1fft(&monomials_jac[0], degree, root, 1);
    for (size_t i = 0; i < degree; ++i) {
        lagrange_jac[i] = g1::group_exponentiation(lagrange_jac[i], domain.domain_inverse);
        g1::__jacobian_to_affine(lagrange_jac[i], lagrange_base_affine[i]);
    }

    /*
        for (size_t i = 0; i < domain.log2_size; ++i) {
            domain.lagrange_base[i] = g1::group_exponentiation(monomials[0], n_inv);
        }
    */
}

} // namespace lagrange_base
}