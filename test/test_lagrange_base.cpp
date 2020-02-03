#include <gtest/gtest.h>

#include <barretenberg/curves/bn254/fq12.hpp>
#include <barretenberg/curves/bn254/g1.hpp>
#include <barretenberg/curves/bn254/g2.hpp>
#include <barretenberg/curves/bn254/scalar_multiplication/scalar_multiplication.hpp>
#include <barretenberg/curves/bn254/pairing.hpp>
#include <barretenberg/io/io.hpp>
#include <barretenberg/polynomials/polynomial.hpp>
#include <barretenberg/types.hpp>

#include <barretenberg/lagrange_base_transformation/lagrange_base.hpp>

using namespace barretenberg;

TEST(lagrange_base, verify_lagrange_base_transformation)
{
    constexpr size_t degree = 64;
    // step 1: create monomial base srs
    // step 2: create lagrange base srs

    // step 3: create polynomial
    // step 4: commit to poly over both reference strings
    // step 5: very correctness
    fr::field_t x = fr::random_element();

    std::array<fr::field_t, degree> monomial_powers;
    monomial_powers[0] = fr::one;
    for (size_t i = 1; i < degree; ++i)
    {
        monomial_powers[i] = fr::mul(monomial_powers[i - 1], x);
    }

    std::array<g1::affine_element, degree * 2> monomial_srs;
    for (size_t i = 0; i < degree; ++i)
    {
        monomial_srs[i] = g1::group_exponentiation(g1::affine_one, monomial_powers[i]);
    }

    std::array<g1::affine_element, degree * 2> lagrange_base_srs;
    lagrange_base::transform_srs(&monomial_srs[0], &lagrange_base_srs[0], degree);

    barretenberg::evaluation_domain domain(degree);
    domain.compute_lookup_table();

    barretenberg::polynomial test_polynomial(degree);
    test_polynomial[0] = fr::one;
    for (size_t i = 1; i < degree; ++i)
    {
        test_polynomial[i] = fr::zero;
    }
    barretenberg::polynomial lagrange_base_polynomial(test_polynomial);

    lagrange_base_polynomial.fft(domain);
    scalar_multiplication::generate_pippenger_point_table(&monomial_srs[0], &monomial_srs[0], degree);
    scalar_multiplication::generate_pippenger_point_table(&lagrange_base_srs[0], &lagrange_base_srs[0], degree);

    g1::element expected = scalar_multiplication::pippenger(&test_polynomial[0], &monomial_srs[0], degree);
    g1::element result = scalar_multiplication::pippenger(&lagrange_base_polynomial[0], &lagrange_base_srs[0], degree);
    expected = g1::normalize(expected);
    result = g1::normalize(result);

    EXPECT_EQ(g1::eq(result, expected), true);
}