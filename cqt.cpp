
#include "cqt.hpp"

#include <math.h>

CArray arrange(size_t N) {
	CArray A(N);
	for (size_t i = 0; i < N; i++) {
		A[i] = 1;
	}
	return A;
}

CArray Cqt::cq_hanning_window(int length)
{
	CArray result(length, 0);

	for (int j = 0; j < length; ++j) 
	{
		float elem = 0.5f * (1 - (cos(j * 2 * PI / (float)(length - 1))));
		result[j] = elem / (float)length;
	}
	return result;
}

Cqt::Cqt(float min_freq, float max_freq, size_t n_bins, float fs) 
{
	_min_freq = min_freq;
	_max_freq = max_freq;
	_n_bins = n_bins;
	_fs = fs;
	_ker = std::vector<CArray>();

	_Q = 1.0 / (pow(2.0, 1.0 / n_bins) - 1.0);
	_fftlen = int(pow(2.0, ceil(log(_Q * fs / min_freq) / log(2.0))));

	int K = ceil(n_bins * log(max_freq / min_freq) / log(2.0));
	double N;
	for (size_t k = K; k > 0; k--)
	{ 
		N = ceil(_Q * fs / (min_freq * pow(2, (k - 1.0) / n_bins)));
		const std::complex<double> j(0.0, 1.0);
		CArray tmpKer = cq_hanning_window(N);
		for (size_t i = 0; i < N; i++)
		{
			tmpKer[i] *= exp(2.0 * PI*j*(double)_Q*(double)i/N)/N;
			fft(tmpKer);
			_ker.push_back(tmpKer);
		}
		std::reverse(_ker.begin(), _ker.end());
	}
			N = ceil(self.Q * fs / (fmin * pow(2, (k - 1.0) / bins)))
			tmpKer = wnd(N) * np.exp(2 * pi * 1j * self.Q * np.arange(N) / N) / N;
		ker = np.fft.fft(tmpKer, self.fftlen)
		# ker = np.select([abs(ker) > self.eps], [ker])
		self.ker += [coo_matrix(ker, dtype = np.complex128)]
		# print 'shape:', hstack(self.ker).tocsc().shape
		self.ker.reverse()
		self.ker = vstack(self.ker).tocsc().transpose().conj() / self.fftlen
}

fftwf_complex* cq_short_time_constq_transform(float* data, int data_length,
	int min_freq, int max_freq, int sample_rate, int bins, int step,
	int* height, int* width) {
	int kernel_height, kernel_width;
	cs_ci* ker = cq_make_kernel(min_freq, max_freq, sample_rate, bins,
		&kernel_height, &kernel_width);
	cs_ci* kern = cs_ci_transpose(ker, 1);

	int max_index = rint(ceil(data_length / (float)kernel_height));

	int indices_size = (max_index - 1) * kernel_height / (double)step + 1;
	int* indices = (int*)malloc(indices_size * sizeof(int));
	for (int j = 0, k = 0; j <= (max_index - 1) * kernel_height; ++k, j += step)
		indices[k] = j;

	fftwf_complex* result = cq_const_q_wrap(data, kern, kernel_height,
		kernel_width, indices, indices_size);

	cs_ci_spfree(ker);
	cs_ci_spfree(kern);
	free(indices);

	*height = kernel_width;
	*width = indices_size;
	return result;
}