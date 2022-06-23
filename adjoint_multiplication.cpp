#include "SHOW.h"
#include "environment.h"

// the stream
constexpr DEG WIDTHIN = 2;
constexpr DEG DEPTHIN = 2;

// the input and output environments
using IN = Environment<WIDTHIN, DEPTHIN>;

typename IN::SHUFFLE_TENSOR shift_down(const typename IN::SHUFFLE_TENSOR& sh, typename IN::SHUFFLE_TENSOR::KEY word)
{
    typename IN::SHUFFLE_TENSOR result(sh), working ;
    while (word.size()) {
        auto letter = word.lparent();
        word = word.rparent();
		for (auto& pr : result) {
            if (pr.key().lparent() == letter)
                working[pr.key().rparent()] = result[pr.key()];
		}
            result.swap(working);
            working.clear();
    }
    return result;
}

// the evaluation of the adjoint operation is worked out below
// <sh,ab>=\sum_{uv=sh}<ua><vb>
//        = <\sum_{uv=sh}<ua>v,b>
// Let T_w(sh) be all the projection of sh onto the part beginning with w with w removed
// \sum_{i} <ki shi,ab> =
//        = \sum_{i} ki<\sum_{uv=shi}<ua>v,b>
//        = \sum_{u} < <ua>T_u(sh), b>
//  The action of the adjoint of tensor multiplication by a multiplication is \sum_u <ua> T_u(sh)     

typename IN::SHUFFLE_TENSOR adjoint_to_multiply(const typename IN::TENSOR& t, typename IN::SHUFFLE_TENSOR sh)
{
	// this implementation is understandable and reliable but repetitive and can be radically accelerated
    IN::SHUFFLE_TENSOR result;
	for (auto& pr : t) {
        result += shift_down(sh, pr.key()) * pr.value();
	}
    return result;
}


int adjoint_multiplication()
{
	auto k = IN::K;
	auto& sbasis = IN::SHUFFLE_TENSOR::basis;
	auto& tbasis = IN::TENSOR::basis;


	IN in;

	std::cout << "Creating a generic shuffle capped to level " << DEPTHIN << "\n\n";

	IN::SHUFFLE_TENSOR sh = in.generic_vector<IN::SHUFFLE_TENSOR>();
	SHOW(sh);

	std::cout << "Creating a generic tensor capped to level " << DEPTHIN << "\n\n";

	IN::TENSOR tensor_ = in.generic_vector<IN::TENSOR>(5000);
    SHOW(tensor_);

	std::cout << "Evaluating the canonical pairing K(shuffle, tensor)" << "\n\n";
	SHOW(IN::K(sh,tensor_));

	std::cout << "Creating the two generic input log signatures \"before\" and \"during\" truncated to level " << DEPTHIN << "\n\n";

	IN::LIE logsig_before = in.generic_vector<IN::LIE>(1000);
	IN::TENSOR sig_before = exp(in.maps_.l2t(logsig_before));
	SHOW(logsig_before);

	IN::LIE logsig_after = in.generic_vector<IN::LIE>(2000);
	IN::TENSOR sig_after = exp(in.maps_.l2t(logsig_after));
	SHOW(logsig_after);

    SHOW(IN::K(sh, sig_before * sig_after));
	SHOW(IN::K(sh, sig_before * sig_after) - IN::K(adjoint_to_multiply(sig_before, sh), sig_after));
	
	return 0;
}