#include "./turbo_composer.hpp"

#include <math.h>

#include "../../curves/bn254/scalar_multiplication/scalar_multiplication.hpp"

#include "../../assert.hpp"
#include "../../curves/bn254/fr.hpp"
#include "../proof_system/widgets/turbo_fixed_base_widget.hpp"
#include "../proof_system/widgets/turbo_arithmetic_widget.hpp"
#include "../proof_system/widgets/turbo_range_widget.hpp"
#include "../proof_system/proving_key/proving_key.hpp"
#include "../proof_system/verification_key/verification_key.hpp"

using namespace barretenberg;

namespace waffle {

TurboComposer::TurboComposer(const size_t size_hint)
    : ComposerBase()
{
    w_l.reserve(size_hint);
    w_r.reserve(size_hint);
    w_o.reserve(size_hint);
    w_4.reserve(size_hint);
    q_m.reserve(size_hint);
    q_1.reserve(size_hint);
    q_2.reserve(size_hint);
    q_3.reserve(size_hint);
    q_4.reserve(size_hint);
    q_arith.reserve(size_hint);
    q_c.reserve(size_hint);
    q_4_next.reserve(size_hint);
    q_ecc_1.reserve(size_hint);
    q_range.reserve(size_hint);

     zero_idx = put_constant_variable(fr::zero);
    // zero_idx = add_variable(barretenberg::fr::zero);
};

void TurboComposer::create_dummy_gate()
{
    gate_flags.push_back(0);
    uint32_t idx = add_variable(fr::to_montgomery_form({{ 1, 0, 0, 0 }}));
    w_l.emplace_back(idx);
    w_r.emplace_back(idx);
    w_o.emplace_back(idx);
    w_4.emplace_back(idx);
    q_arith.emplace_back(fr::zero);
    q_4.emplace_back(fr::zero);
    q_4_next.emplace_back(fr::zero);
    q_ecc_1.emplace_back(fr::zero);
    q_m.emplace_back(fr::zero);
    q_1.emplace_back(fr::zero);
    q_2.emplace_back(fr::zero);
    q_3.emplace_back(fr::zero);
    q_c.emplace_back(fr::zero);
    q_range.emplace_back(fr::zero);
    ++n;
}

void TurboComposer::create_add_gate(const add_triple& in)
{
    gate_flags.push_back(0);
    w_l.emplace_back(in.a);
    w_r.emplace_back(in.b);
    w_o.emplace_back(in.c);
    w_4.emplace_back(zero_idx);
    q_m.emplace_back(fr::zero);
    q_1.emplace_back(in.a_scaling);
    q_2.emplace_back(in.b_scaling);
    q_3.emplace_back(in.c_scaling);
    q_c.emplace_back(in.const_scaling);
    q_arith.emplace_back(fr::one);
    q_4.emplace_back(fr::zero);
    q_4_next.emplace_back(fr::zero);
    q_ecc_1.emplace_back(fr::zero);
    q_range.emplace_back(fr::zero);

    epicycle left{ static_cast<uint32_t>(n), WireType::LEFT };
    epicycle right{ static_cast<uint32_t>(n), WireType::RIGHT };
    epicycle out{ static_cast<uint32_t>(n), WireType::OUTPUT };
    epicycle fourth{ static_cast<uint32_t>(n), WireType::FOURTH };

    ASSERT(wire_epicycles.size() > in.a);
    ASSERT(wire_epicycles.size() > in.b);
    ASSERT(wire_epicycles.size() > in.c);
    ASSERT(wire_epicycles.size() > zero_idx);

    wire_epicycles[static_cast<size_t>(in.a)].emplace_back(left);
    wire_epicycles[static_cast<size_t>(in.b)].emplace_back(right);
    wire_epicycles[static_cast<size_t>(in.c)].emplace_back(out);
    wire_epicycles[static_cast<size_t>(zero_idx)].emplace_back(fourth);

    ++n;
}

void TurboComposer::create_big_add_gate(const add_quad& in)
{
    gate_flags.push_back(0);
    w_l.emplace_back(in.a);
    w_r.emplace_back(in.b);
    w_o.emplace_back(in.c);
    w_4.emplace_back(in.d);
    q_m.emplace_back(fr::zero);
    q_1.emplace_back(in.a_scaling);
    q_2.emplace_back(in.b_scaling);
    q_3.emplace_back(in.c_scaling);
    q_c.emplace_back(in.const_scaling);
    q_arith.emplace_back(fr::one);
    q_4.emplace_back(in.d_scaling);
    q_4_next.emplace_back(fr::zero);
    q_ecc_1.emplace_back(fr::zero);
    q_range.emplace_back(fr::zero);

    epicycle left{ static_cast<uint32_t>(n), WireType::LEFT };
    epicycle right{ static_cast<uint32_t>(n), WireType::RIGHT };
    epicycle out{ static_cast<uint32_t>(n), WireType::OUTPUT };
    epicycle fourth{ static_cast<uint32_t>(n), WireType::FOURTH };

    ASSERT(wire_epicycles.size() > in.a);
    ASSERT(wire_epicycles.size() > in.b);
    ASSERT(wire_epicycles.size() > in.c);
    ASSERT(wire_epicycles.size() > in.d);

    wire_epicycles[static_cast<size_t>(in.a)].emplace_back(left);
    wire_epicycles[static_cast<size_t>(in.b)].emplace_back(right);
    wire_epicycles[static_cast<size_t>(in.c)].emplace_back(out);
    wire_epicycles[static_cast<size_t>(in.d)].emplace_back(fourth);

    ++n;
}

void TurboComposer::create_mul_gate(const mul_triple& in)
{
    gate_flags.push_back(0);
    add_gate_flag(gate_flags.size() - 1, GateFlags::FIXED_LEFT_WIRE);
    add_gate_flag(gate_flags.size() - 1, GateFlags::FIXED_RIGHT_WIRE);
    w_l.emplace_back(in.a);
    w_r.emplace_back(in.b);
    w_o.emplace_back(in.c);
    w_4.emplace_back(zero_idx);
    q_m.emplace_back(in.mul_scaling);
    q_1.emplace_back(fr::zero);
    q_2.emplace_back(fr::zero);
    q_3.emplace_back(in.c_scaling);
    q_c.emplace_back(in.const_scaling);
    q_arith.emplace_back(fr::one);
    q_4.emplace_back(fr::zero);
    q_4_next.emplace_back(fr::zero);
    q_ecc_1.emplace_back(fr::zero);
    q_range.emplace_back(fr::zero);

    epicycle left{ static_cast<uint32_t>(n), WireType::LEFT };
    epicycle right{ static_cast<uint32_t>(n), WireType::RIGHT };
    epicycle out{ static_cast<uint32_t>(n), WireType::OUTPUT };

    ASSERT(wire_epicycles.size() > in.a);
    ASSERT(wire_epicycles.size() > in.b);
    ASSERT(wire_epicycles.size() > in.c);
    ASSERT(wire_epicycles.size() > zero_idx);

    wire_epicycles[static_cast<size_t>(in.a)].emplace_back(left);
    wire_epicycles[static_cast<size_t>(in.b)].emplace_back(right);
    wire_epicycles[static_cast<size_t>(in.c)].emplace_back(out);

    ++n;
}

void TurboComposer::create_bool_gate(const uint32_t variable_index)
{
    gate_flags.push_back(0);
    add_gate_flag(gate_flags.size() - 1, GateFlags::FIXED_LEFT_WIRE);
    add_gate_flag(gate_flags.size() - 1, GateFlags::FIXED_RIGHT_WIRE);
    w_l.emplace_back(variable_index);
    w_r.emplace_back(variable_index);
    w_o.emplace_back(variable_index);
    w_4.emplace_back(zero_idx);
    q_arith.emplace_back(fr::one);
    q_4.emplace_back(fr::zero);
    q_4_next.emplace_back(fr::zero);
    q_ecc_1.emplace_back(fr::zero);
    q_range.emplace_back(fr::zero);

    q_m.emplace_back(fr::one);
    q_1.emplace_back(fr::zero);
    q_2.emplace_back(fr::zero);
    q_3.emplace_back(fr::neg_one());
    q_c.emplace_back(fr::zero);

    epicycle left{ static_cast<uint32_t>(n), WireType::LEFT };
    epicycle right{ static_cast<uint32_t>(n), WireType::RIGHT };
    epicycle out{ static_cast<uint32_t>(n), WireType::OUTPUT };

    ASSERT(wire_epicycles.size() > variable_index);
    wire_epicycles[static_cast<size_t>(variable_index)].emplace_back(left);
    wire_epicycles[static_cast<size_t>(variable_index)].emplace_back(right);
    wire_epicycles[static_cast<size_t>(variable_index)].emplace_back(out);

    ++n;
}

void TurboComposer::create_poly_gate(const poly_triple& in)
{
    gate_flags.push_back(0);
    add_gate_flag(gate_flags.size() - 1, GateFlags::FIXED_LEFT_WIRE);
    add_gate_flag(gate_flags.size() - 1, GateFlags::FIXED_RIGHT_WIRE);
    w_l.emplace_back(in.a);
    w_r.emplace_back(in.b);
    w_o.emplace_back(in.c);
    w_4.emplace_back(zero_idx);
    q_m.emplace_back(in.q_m);
    q_1.emplace_back(in.q_l);
    q_2.emplace_back(in.q_r);
    q_3.emplace_back(in.q_o);
    q_c.emplace_back(in.q_c);
    q_range.emplace_back(fr::zero);

    q_arith.emplace_back(fr::one);
    q_4.emplace_back(fr::zero);
    q_4_next.emplace_back(fr::zero);
    q_ecc_1.emplace_back(fr::zero);

    epicycle left{ static_cast<uint32_t>(n), WireType::LEFT };
    epicycle right{ static_cast<uint32_t>(n), WireType::RIGHT };
    epicycle out{ static_cast<uint32_t>(n), WireType::OUTPUT };

    ASSERT(wire_epicycles.size() > in.a);
    ASSERT(wire_epicycles.size() > in.b);
    ASSERT(wire_epicycles.size() > in.c);
    ASSERT(wire_epicycles.size() > zero_idx);

    wire_epicycles[static_cast<size_t>(in.a)].emplace_back(left);
    wire_epicycles[static_cast<size_t>(in.b)].emplace_back(right);
    wire_epicycles[static_cast<size_t>(in.c)].emplace_back(out);

    ++n;
}

void TurboComposer::create_fixed_group_add_gate(const fixed_group_add_quad& in)
{
    gate_flags.push_back(0);
    w_l.emplace_back(in.a);
    w_r.emplace_back(in.b);
    w_o.emplace_back(in.c);
    w_4.emplace_back(in.d);

    q_arith.emplace_back(fr::zero);
    q_4.emplace_back(fr::zero);
    q_4_next.emplace_back(fr::zero);
    q_m.emplace_back(fr::zero);
    q_c.emplace_back(fr::zero);
    q_range.emplace_back(fr::zero);

    q_1.emplace_back(in.q_x_1);
    q_2.emplace_back(in.q_x_2);
    q_3.emplace_back(in.q_y_1);
    q_ecc_1.emplace_back(in.q_y_2);

    epicycle left{ static_cast<uint32_t>(n), WireType::LEFT };
    epicycle right{ static_cast<uint32_t>(n), WireType::RIGHT };
    epicycle out{ static_cast<uint32_t>(n), WireType::OUTPUT };
    epicycle fourth{ static_cast<uint32_t>(n), WireType::FOURTH };

    ASSERT(wire_epicycles.size() > in.a);
    ASSERT(wire_epicycles.size() > in.b);
    ASSERT(wire_epicycles.size() > in.c);
    ASSERT(wire_epicycles.size() > in.d);

    wire_epicycles[static_cast<size_t>(in.a)].emplace_back(left);
    wire_epicycles[static_cast<size_t>(in.b)].emplace_back(right);
    wire_epicycles[static_cast<size_t>(in.c)].emplace_back(out);
    wire_epicycles[static_cast<size_t>(in.d)].emplace_back(fourth);

    ++n;
}

void TurboComposer::create_fixed_group_add_gate_with_init(const fixed_group_add_quad& in, const fixed_group_init_quad& init)
{
    gate_flags.push_back(0);
    w_l.emplace_back(in.a);
    w_r.emplace_back(in.b);
    w_o.emplace_back(in.c);
    w_4.emplace_back(in.d);

    q_arith.emplace_back(fr::zero);
    q_4.emplace_back(init.q_x_1);
    q_4_next.emplace_back(init.q_x_2);
    q_m.emplace_back(init.q_y_1);
    q_c.emplace_back(init.q_y_2);
    q_range.emplace_back(fr::zero);

    q_1.emplace_back(in.q_x_1);
    q_2.emplace_back(in.q_x_2);
    q_3.emplace_back(in.q_y_1);
    q_ecc_1.emplace_back(in.q_y_2);

    epicycle left{ static_cast<uint32_t>(n), WireType::LEFT };
    epicycle right{ static_cast<uint32_t>(n), WireType::RIGHT };
    epicycle out{ static_cast<uint32_t>(n), WireType::OUTPUT };
    epicycle fourth{ static_cast<uint32_t>(n), WireType::FOURTH };

    ASSERT(wire_epicycles.size() > in.a);
    ASSERT(wire_epicycles.size() > in.b);
    ASSERT(wire_epicycles.size() > in.c);
    ASSERT(wire_epicycles.size() > in.d);

    wire_epicycles[static_cast<size_t>(in.a)].emplace_back(left);
    wire_epicycles[static_cast<size_t>(in.b)].emplace_back(right);
    wire_epicycles[static_cast<size_t>(in.c)].emplace_back(out);
    wire_epicycles[static_cast<size_t>(in.d)].emplace_back(fourth);

    ++n;  
}

void TurboComposer::fix_witness(const uint32_t witness_index, const barretenberg::fr::field_t& witness_value)
{
    gate_flags.push_back(0);

    w_l.emplace_back(witness_index);
    w_r.emplace_back(zero_idx);
    w_o.emplace_back(zero_idx);
    w_4.emplace_back(zero_idx);
    q_m.emplace_back(fr::zero);
    q_1.emplace_back(fr::one);
    q_2.emplace_back(fr::zero);
    q_3.emplace_back(fr::zero);
    q_c.emplace_back(fr::neg(witness_value));
    q_arith.emplace_back(fr::one);
    q_4.emplace_back(fr::zero);
    q_4_next.emplace_back(fr::zero);
    q_ecc_1.emplace_back(fr::zero);
    q_range.emplace_back(fr::zero);

    epicycle left{ static_cast<uint32_t>(n), WireType::LEFT };

    ASSERT(wire_epicycles.size() > witness_index);
    ASSERT(wire_epicycles.size() > zero_idx);
    ASSERT(wire_epicycles.size() > zero_idx);
    wire_epicycles[static_cast<size_t>(witness_index)].emplace_back(left);

    ++n;
}

std::vector<uint32_t> TurboComposer::create_range_constraint(const uint32_t witness_index, const size_t num_bits)
{
    ASSERT(static_cast<uint32_t>(variables.size()) > witness_index);
    ASSERT(((num_bits >> 1U) << 1U) == num_bits);

    /*
     * The range constraint accumulates base 4 values into a sum.
     * We do this by evaluating a kind of 'raster scan', where we compare adjacent elements
     * and validate that their differences map to a base for value  * 
     * Let's say that we want to perform a 32-bit range constraint in 'x'.
     * We can represent x via 16 constituent base-4 'quads' {q_0, ..., q_15}:
     * 
     *      15           
     *      ===          
     *      \          i
     * x =  /    q  . 4 
     *      ===   i       
     *     i = 0
     * 
     * In program memory, we place an accumulating base-4 sum of x {a_0, ..., a_15}, where
     * 
     *         i                  
     *        ===                 
     *        \                  j
     * a   =  /    q         .  4 
     *  i     ===   (15 - j)              
     *       j = 0   
     *
     *
     * From this, we can use our range transition constraint to validate that
     * 
     * 
     *  a      - 4 . a  ϵ [0, 1, 2, 3]
     *   i + 1        i
     * 
     * 
     * We place our accumulating sums in program memory in the following sequence: 
     * 
     * +-----+-----+-----+-----+
     * |  A  |  B  |  C  |  D  |
     * +-----+-----+-----+-----+
     * | a3  | a2  | a1  | 0   |
     * | a7  | a6  | a5  | a4  |
     * | a11 | a10 | a9  | a8  |
     * | a15 | a14 | a13 | a12 |
     * | --- | --- | --- | a16 |
     * +-----+-----+-----+-----+
     *
     * Our range transition constraint on row 'i'
     * performs our base-4 range check on the follwing pairs:
     * 
     * (D_{i}, C_{i}), (C_{i}, B_{i}), (B_{i}, A_{i}), (A_{i}, D_{i+1})
     * 
     * We need to start our raster scan at zero, so we simplify matters and just force the first value
     * to be zero.
     * 
     * The output will be in the 4th column of an otherwise unused row. Assuming this row can
     * be used for a width-3 standard gate, the total number of gates for an n-bit range constraint
     * is (n / 8) gates
     *
     **/
 
    const fr::field_t witness_value = fr::from_montgomery_form(variables[witness_index]);

    // one gate accmulates 4 quads, or 8 bits.
    // # gates = (bits / 8)
    size_t num_quad_gates = (num_bits >> 3);
    
    num_quad_gates = (num_quad_gates << 3 == num_bits) ? num_quad_gates : num_quad_gates + 1;

    // hmm
    std::vector<uint32_t>* wires[4]{ &w_4, &w_o, &w_r, &w_l };

    // hmmm
    WireType wire_types[4]{ WireType::FOURTH, WireType::OUTPUT, WireType::RIGHT, WireType::LEFT };

    const size_t num_quads = (num_quad_gates << 2);
    const size_t forced_zero_threshold = 1 + (((num_quads << 1) - num_bits) >> 1);
    std::vector<uint32_t> accumulators;
    fr::field_t accumulator = fr::zero;

    for (size_t i = 0; i < num_quads + 1; ++i)
    {
        const size_t gate_index = n + (i / 4);
        uint32_t accumulator_index;
        if (i < forced_zero_threshold)
        {
            accumulator_index = zero_idx;
        }
        else
        {
            const size_t bit_index = (num_quads - i) << 1;
            const uint64_t quad = static_cast<uint64_t>(fr::get_bit(witness_value, bit_index)) + 2ULL * static_cast<uint64_t>(fr::get_bit(witness_value, bit_index + 1));
            const fr::field_t quad_element = fr::to_montgomery_form({ quad, 0, 0, 0 });
            fr::__add(accumulator, accumulator, accumulator);
            fr::__add(accumulator, accumulator, accumulator);
            fr::__add(accumulator, quad_element, accumulator);

            accumulator_index = add_variable(accumulator);
            accumulators.emplace_back(accumulator_index);
        }

        // hmmmm
        (*(wires + (i & 3)))->emplace_back(accumulator_index);
        const size_t wire_index = i & 3;

        wire_epicycles[accumulator_index].emplace_back(epicycle(static_cast<uint32_t>(gate_index), wire_types[wire_index]));
    }
    size_t used_gates = (num_quads + 1) / 4;

    // TODO: handle partially used gates. For now just set them to be zero
    if (used_gates * 4 != (num_quads + 1))
    {
        ++used_gates;
    }

    for (size_t i = 0; i < used_gates; ++i)
    {
        q_m.emplace_back(fr::zero);
        q_1.emplace_back(fr::zero);
        q_2.emplace_back(fr::zero);
        q_3.emplace_back(fr::zero);
        q_c.emplace_back(fr::zero);
        q_arith.emplace_back(fr::zero);
        q_4.emplace_back(fr::zero);
        q_4_next.emplace_back(fr::zero);
        q_ecc_1.emplace_back(fr::zero);

        q_range.emplace_back(fr::one);
    }

    q_range[q_range.size() - 1] = fr::zero;

    w_l.emplace_back(zero_idx);
    w_r.emplace_back(zero_idx);
    w_o.emplace_back(zero_idx);

    assert_equal(accumulators[accumulators.size() - 1], witness_index);
    accumulators[accumulators.size() - 1] = witness_index;

    n += used_gates;
    return accumulators;
}

uint32_t TurboComposer::put_constant_variable(const barretenberg::fr::field_t& variable)
{
    if (constant_variables.count(variable) == 1)
    {
        return constant_variables.at(variable);
    }
    else
    {
        uint32_t variable_index = add_variable(variable);
        fix_witness(variable_index, variable);
        constant_variables.insert({ variable, variable_index });
        return variable_index;
    }
}

std::shared_ptr<proving_key> TurboComposer::compute_proving_key()
{
    if (computed_proving_key)
    {
        return circuit_proving_key;
    }
    create_dummy_gate();
    ASSERT(wire_epicycles.size() == variables.size());
    ASSERT(n == q_m.size());
    ASSERT(n == q_1.size());
    ASSERT(n == q_2.size());
    ASSERT(n == q_3.size());
    ASSERT(n == q_3.size());
    ASSERT(n == q_4.size());
    ASSERT(n == q_4_next.size());
    ASSERT(n == q_arith.size());
    ASSERT(n == q_ecc_1.size());
    ASSERT(n == q_range.size());

    const size_t total_num_gates = n + public_inputs.size();

    size_t log2_n = static_cast<size_t>(log2(total_num_gates + 1));
    if ((1UL << log2_n) != (total_num_gates + 1)) {
        ++log2_n;
    }
    size_t new_n = 1UL << log2_n;

    for (size_t i = total_num_gates; i < new_n; ++i) {
        q_m.emplace_back(fr::zero);
        q_1.emplace_back(fr::zero);
        q_2.emplace_back(fr::zero);
        q_3.emplace_back(fr::zero);
        q_c.emplace_back(fr::zero);
        q_4.emplace_back(fr::zero);
        q_4_next.emplace_back(fr::zero);
        q_arith.emplace_back(fr::zero);
        q_ecc_1.emplace_back(fr::zero);
        q_range.emplace_back(fr::zero);
    }

    for (size_t i = 0; i < public_inputs.size(); ++i)
    {
        epicycle left{ static_cast<uint32_t>(i - public_inputs.size()), WireType::LEFT };
        wire_epicycles[static_cast<size_t>(public_inputs[i])].emplace_back(left);
    }
    circuit_proving_key = std::make_shared<proving_key>(new_n, public_inputs.size());
    
    polynomial poly_q_m(new_n);
    polynomial poly_q_c(new_n);
    polynomial poly_q_1(new_n);
    polynomial poly_q_2(new_n);
    polynomial poly_q_3(new_n);
    polynomial poly_q_4(new_n);
    polynomial poly_q_4_next(new_n);
    polynomial poly_q_arith(new_n);
    polynomial poly_q_ecc_1(new_n);
    polynomial poly_q_range(new_n);


    for (size_t i = 0; i < public_inputs.size(); ++i)
    {
        poly_q_m[i] = fr::zero;
        poly_q_1[i] = fr::zero;
        poly_q_2[i] = fr::zero;
        poly_q_3[i] = fr::zero;
        poly_q_4[i] = fr::zero;
        poly_q_4_next[i] = fr::zero;
        poly_q_arith[i] = fr::zero;
        poly_q_ecc_1[i] = fr::zero;
        poly_q_c[i] = fr::zero;
        poly_q_range[i] = fr::zero;
    }

    for (size_t i = public_inputs.size(); i < new_n; ++i)
    {
        poly_q_m[i] = q_m[i - public_inputs.size()];
        poly_q_1[i] = q_1[i - public_inputs.size()];
        poly_q_2[i] = q_2[i - public_inputs.size()];
        poly_q_3[i] = q_3[i - public_inputs.size()];
        poly_q_c[i] = q_c[i - public_inputs.size()];
        poly_q_4[i] = q_4[i - public_inputs.size()];
        poly_q_4_next[i] = q_4_next[i - public_inputs.size()];
        poly_q_arith[i] = q_arith[i - public_inputs.size()];
        poly_q_ecc_1[i] = q_ecc_1[i - public_inputs.size()];
        poly_q_range[i] = q_range[i - public_inputs.size()];
    }

    poly_q_1.ifft(circuit_proving_key->small_domain);
    poly_q_2.ifft(circuit_proving_key->small_domain);
    poly_q_3.ifft(circuit_proving_key->small_domain);
    poly_q_4.ifft(circuit_proving_key->small_domain);
    poly_q_4_next.ifft(circuit_proving_key->small_domain);
    poly_q_m.ifft(circuit_proving_key->small_domain);
    poly_q_c.ifft(circuit_proving_key->small_domain);
    poly_q_arith.ifft(circuit_proving_key->small_domain);
    poly_q_ecc_1.ifft(circuit_proving_key->small_domain);
    poly_q_range.ifft(circuit_proving_key->small_domain);

    polynomial poly_q_1_fft(poly_q_1, new_n * 4);
    polynomial poly_q_2_fft(poly_q_2, new_n * 4);
    polynomial poly_q_3_fft(poly_q_3, new_n * 4);
    polynomial poly_q_4_fft(poly_q_4, new_n * 4);
    polynomial poly_q_4_next_fft(poly_q_4_next, new_n * 4);
    polynomial poly_q_m_fft(poly_q_m, new_n * 4);
    polynomial poly_q_c_fft(poly_q_c, new_n * 4);
    polynomial poly_q_arith_fft(poly_q_arith, new_n * 4);
    polynomial poly_q_ecc_1_fft(poly_q_ecc_1, new_n * 4);
    polynomial poly_q_range_fft(poly_q_range, new_n * 4);

    poly_q_1_fft.coset_fft(circuit_proving_key->large_domain);
    poly_q_2_fft.coset_fft(circuit_proving_key->large_domain);
    poly_q_3_fft.coset_fft(circuit_proving_key->large_domain);
    poly_q_4_fft.coset_fft(circuit_proving_key->large_domain);
    poly_q_4_next_fft.coset_fft(circuit_proving_key->large_domain);
    poly_q_m_fft.coset_fft(circuit_proving_key->large_domain);
    poly_q_c_fft.coset_fft(circuit_proving_key->large_domain);
    poly_q_arith_fft.coset_fft(circuit_proving_key->large_domain);
    poly_q_ecc_1_fft.coset_fft(circuit_proving_key->large_domain);
    poly_q_range_fft.coset_fft(circuit_proving_key->large_domain);

    circuit_proving_key->constraint_selectors.insert({ "q_m", std::move(poly_q_m )});
    circuit_proving_key->constraint_selectors.insert({ "q_c", std::move(poly_q_c )});
    circuit_proving_key->constraint_selectors.insert({ "q_arith", std::move(poly_q_arith )});
    circuit_proving_key->constraint_selectors.insert({ "q_ecc_1", std::move(poly_q_ecc_1 )});
    circuit_proving_key->constraint_selectors.insert({ "q_1", std::move(poly_q_1 )});
    circuit_proving_key->constraint_selectors.insert({ "q_2", std::move(poly_q_2 )});
    circuit_proving_key->constraint_selectors.insert({ "q_3", std::move(poly_q_3 )});
    circuit_proving_key->constraint_selectors.insert({ "q_4", std::move(poly_q_4 )});
    circuit_proving_key->constraint_selectors.insert({ "q_4_next", std::move(poly_q_4_next )});
    circuit_proving_key->constraint_selectors.insert({ "q_range", std::move(poly_q_range )});

    circuit_proving_key->constraint_selector_ffts.insert({ "q_m_fft", std::move(poly_q_m_fft )});
    circuit_proving_key->constraint_selector_ffts.insert({ "q_c_fft", std::move(poly_q_c_fft )});
    circuit_proving_key->constraint_selector_ffts.insert({ "q_arith_fft", std::move(poly_q_arith_fft )});
    circuit_proving_key->constraint_selector_ffts.insert({ "q_ecc_1_fft", std::move(poly_q_ecc_1_fft )});
    circuit_proving_key->constraint_selector_ffts.insert({ "q_1_fft", std::move(poly_q_1_fft )});
    circuit_proving_key->constraint_selector_ffts.insert({ "q_2_fft", std::move(poly_q_2_fft )});
    circuit_proving_key->constraint_selector_ffts.insert({ "q_3_fft", std::move(poly_q_3_fft )});
    circuit_proving_key->constraint_selector_ffts.insert({ "q_4_fft", std::move(poly_q_4_fft )});
    circuit_proving_key->constraint_selector_ffts.insert({ "q_4_next_fft", std::move(poly_q_4_next_fft )});
    circuit_proving_key->constraint_selector_ffts.insert({ "q_range_fft", std::move(poly_q_range_fft )});

    compute_sigma_permutations<4>(circuit_proving_key.get());
    computed_proving_key = true;
    return circuit_proving_key;
}

std::shared_ptr<verification_key> TurboComposer::compute_verification_key()
{
    if (computed_verification_key)
    {
        return circuit_verification_key;
    }
    if (!computed_proving_key)
    {
        compute_proving_key();
    }

    std::array<fr::field_t*, 13> poly_coefficients;
    poly_coefficients[0] = circuit_proving_key->constraint_selectors.at("q_1").get_coefficients();
    poly_coefficients[1] = circuit_proving_key->constraint_selectors.at("q_2").get_coefficients();
    poly_coefficients[2] = circuit_proving_key->constraint_selectors.at("q_3").get_coefficients();
    poly_coefficients[3] = circuit_proving_key->constraint_selectors.at("q_4").get_coefficients();
    poly_coefficients[4] = circuit_proving_key->constraint_selectors.at("q_4_next").get_coefficients();
    poly_coefficients[5] = circuit_proving_key->constraint_selectors.at("q_m").get_coefficients();
    poly_coefficients[6] = circuit_proving_key->constraint_selectors.at("q_c").get_coefficients();
    poly_coefficients[7] = circuit_proving_key->constraint_selectors.at("q_arith").get_coefficients();
    poly_coefficients[8] = circuit_proving_key->constraint_selectors.at("q_ecc_1").get_coefficients();
    poly_coefficients[9] = circuit_proving_key->constraint_selectors.at("q_range").get_coefficients();

    poly_coefficients[10] = circuit_proving_key->permutation_selectors.at("sigma_1").get_coefficients();
    poly_coefficients[11] = circuit_proving_key->permutation_selectors.at("sigma_2").get_coefficients();
    poly_coefficients[12] = circuit_proving_key->permutation_selectors.at("sigma_3").get_coefficients();

    std::vector<barretenberg::g1::affine_element> commitments;
    commitments.resize(13);

    for (size_t i = 0; i < 13; ++i) {
        g1::jacobian_to_affine(
            scalar_multiplication::pippenger(poly_coefficients[i], circuit_proving_key->reference_string.monomials, circuit_proving_key->n),
            commitments[i]);
    }

    circuit_verification_key = std::make_shared<verification_key>(circuit_proving_key->n);

    circuit_verification_key->constraint_selectors.insert({ "Q_1", commitments[0] });
    circuit_verification_key->constraint_selectors.insert({ "Q_2", commitments[1] });
    circuit_verification_key->constraint_selectors.insert({ "Q_3", commitments[2] });
    circuit_verification_key->constraint_selectors.insert({ "Q_4", commitments[3] });
    circuit_verification_key->constraint_selectors.insert({ "Q_4_NEXT", commitments[4] });
    circuit_verification_key->constraint_selectors.insert({ "Q_M", commitments[5] });
    circuit_verification_key->constraint_selectors.insert({ "Q_C", commitments[6] });
    circuit_verification_key->constraint_selectors.insert({ "Q_ARITH", commitments[7] });
    circuit_verification_key->constraint_selectors.insert({ "Q_ECC_1", commitments[8] });
    circuit_verification_key->constraint_selectors.insert({ "Q_RANGE", commitments[9] });

    circuit_verification_key->permutation_selectors.insert({ "SIGMA_1", commitments[10] });
    circuit_verification_key->permutation_selectors.insert({ "SIGMA_2", commitments[11] });
    circuit_verification_key->permutation_selectors.insert({ "SIGMA_3", commitments[12] });

    computed_verification_key = true;
    return circuit_verification_key;
}

std::shared_ptr<program_witness> TurboComposer::compute_witness()
{
    if (computed_witness)
    {
        return witness;
    }
    const size_t total_num_gates = n + public_inputs.size();
    size_t log2_n = static_cast<size_t>(log2(total_num_gates + 1));
    if ((1UL << log2_n) != (total_num_gates + 1)) {
        ++log2_n;
    }
    size_t new_n = 1UL << log2_n;

    for (size_t i = total_num_gates; i < new_n; ++i) {
        w_l.emplace_back(zero_idx);
        w_r.emplace_back(zero_idx);
        w_o.emplace_back(zero_idx);
        w_4.emplace_back(zero_idx);
    }

    polynomial poly_w_1(new_n);
    polynomial poly_w_2(new_n);
    polynomial poly_w_3(new_n);
    polynomial poly_w_4(new_n);

    for (size_t i = 0; i < public_inputs.size(); ++i)
    {
        fr::__copy(variables[public_inputs[i]], poly_w_1[i]);
        fr::__copy(fr::zero, poly_w_2[i]);
        fr::__copy(fr::zero, poly_w_3[i]);
        fr::__copy(fr::zero, poly_w_4[i]);
    }
    for (size_t i = public_inputs.size(); i < new_n; ++i) {
        fr::__copy(variables[w_l[i - public_inputs.size()]], poly_w_1.at(i));
        fr::__copy(variables[w_r[i - public_inputs.size()]], poly_w_2.at(i));
        fr::__copy(variables[w_o[i - public_inputs.size()]], poly_w_3.at(i));
        fr::__copy(variables[w_4[i - public_inputs.size()]], poly_w_4.at(i));
    }

    witness = std::make_shared<program_witness>();
    witness->wires.insert({ "w_1", std::move(poly_w_1) });
    witness->wires.insert({ "w_2", std::move(poly_w_2) });
    witness->wires.insert({ "w_3", std::move(poly_w_3) });
    witness->wires.insert({ "w_4", std::move(poly_w_4) });

    computed_witness = true;
    return witness;
}

TurboProver TurboComposer::preprocess()
{
    compute_proving_key();
    compute_witness();

    TurboProver output_state(circuit_proving_key, witness, create_manifest(public_inputs.size()));

    std::unique_ptr<ProverTurboFixedBaseWidget> fixed_base_widget = std::make_unique<ProverTurboFixedBaseWidget>(circuit_proving_key.get(), witness.get());
    std::unique_ptr<ProverTurboRangeWidget> range_widget = std::make_unique<ProverTurboRangeWidget>(circuit_proving_key.get(), witness.get());

    output_state.widgets.emplace_back(std::move(fixed_base_widget));
    output_state.widgets.emplace_back(std::move(range_widget));

    return output_state;
}
} // namespace waffle