#include "SHOW.h"
#include "environment.h"

// the stream
constexpr DEG WIDTHIN = 2;
constexpr DEG DEPTHIN = 2;

// the input and output environments
using IN = Environment<WIDTHIN, DEPTHIN>;

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
	SHOW(IN::K(sh, sig_before * sig_after) - IN::K(IN::adjoint_to_multiply(sig_before, sh), sig_after));

	IN::TENSOR sig(sig_before * sig_after);
    SHOW(IN::K(sh, antipode(sig_after)) - IN::K(IN::adjoint_to_multiply(antipode(sig),sh), sig_before));

	return 0;
}