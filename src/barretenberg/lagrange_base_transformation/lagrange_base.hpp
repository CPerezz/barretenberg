#pragma once 

#include "../curves/bn254/g1.hpp"
#include "../polynomials/evaluation_domain.hpp"
/*
#include <barretenberg/curves/bn254/scalar_multiplication/scalar_multiplication.hpp>
#include <barretenberg/curves/bn254/pairing.hpp>
#include <barretenberg/io/io.hpp>
#include <barretenberg/polynomials/polynomial.hpp>
#include <barretenberg/types.hpp>
*/

namespace barretenberg
{
namespace lagrange_base
{

    void transform_srs(g1::affine_element*, g1::affine_element*, const size_t);
}
}