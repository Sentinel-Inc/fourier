//
// Created by M.Kojro on 1/6/2022.
//

#ifndef FOURIER_TRANSFORM_SRC_FOURIER_H_
#define FOURIER_TRANSFORM_SRC_FOURIER_H_
#include <cassert>
#include <complex>
#include <vector>
namespace fourier {

///
/// \param signal
/// \return
std::vector<std::complex<double>> dft(const std::vector<double>& signal);

///
/// \param signal
/// \return
std::vector<double> ift(const std::vector<std::complex<double>>& signal);

#if FOURIER_WITH_FFT

///
/// \param signal
/// \return
double fft(const std::vector<std::complex<double>>& signal);

#endif
		
///
/// \param signal
/// \param x
/// \param y
/// \return
std::complex<double> dft2d(const std::vector<std::vector<double>> &signal,
                           int x, int y);

/// \note  not yet tested
/// \param signal
/// \return
std::vector<std::vector<double>>
dft2d(const std::vector<std::vector<double>> &signal);

#ifdef FOURIER_IMPLEMENTATIONS

static const double kPi = 3.1415;

std::complex<double> dft_impl(const std::vector<double>& signal, double freq) {
  using namespace std::complex_literals;

  std::complex<double> result = 0;
  double N = signal.size();
  for (int i = 0; i < signal.size(); i++) {
	double n = i;
	double k = freq;
	double x = signal[i];
	result += x * exp(-(1i * 2. * kPi * n * k) / N);
  }
  return 1. / N * result;
}

std::vector<std::complex<double>> dft(const std::vector<double>& signal) {
  std::vector<std::complex<double>> tramsformed(signal.size());

  for (int x = 0; x < signal.size(); x++) {
	tramsformed[x] = dft_impl(signal, x);
  }
  return tramsformed;
}

double ift_impl(const std::vector<std::complex<double>>& signal, double x) {
  using namespace std::complex_literals;

  std::complex<double> result = 0;
  for (int i = 0; i < signal.size(); i++) {
	double k = i;
	double N = signal.size();
	std::complex<double> fk = signal[i];
	result += exp((2.0i * kPi * k * x) / N) * fk;
  }

  // The imaginary component is equal to 0
  return result.real();
}

std::vector<double> ift(const std::vector<std::complex<double>>& signal) {
  std::vector<double> reversed(signal.size());
  for (int x = 0; x < signal.size(); x++) { reversed[x] = ift_impl(signal, x); }
  return reversed;
}

#if FOURIER_WITH_FFT

template<class T>
std::vector<T> concat(const std::vector<std::vector<T>>& vecs) {
  std::vector<T> result;

  size_t result_sz = 0;
  for (auto& v : vecs) { result_sz += v.size(); }

  result.reserve(result_sz);
  for (auto& v : vecs) {
	for (auto& elem : v) { result.push_back(elem); }
  }

  return result;
}

std::vector<std::complex<double>> mul(
	const std::vector<std::complex<double>>& a,
	const std::vector<std::complex<double>>& b) {

  assert(a.size() == b.size());
  std::vector<std::complex<double>> result(a.size());
  for (size_t i = 0; i < a.size(); i++) { result[i] = a[i] * b[i]; }
  return result;
}

std::vector<std::complex<double>> add(
	const std::vector<std::complex<double>>& a,
	const std::vector<std::complex<double>>& b) {

  assert(a.size() == b.size());
  std::vector<std::complex<double>> result(a.size());
  for (size_t i = 0; i < a.size(); i++) { result[i] = a[i] + b[i]; }
  return result;
}

std::vector<std::complex<double>> fft(const std::vector<double>& signal) {
  using namespace std::complex_literals;

  size_t N = signal.size();
  size_t cutoff_size = 16;

  assert(N % 2 == 0 && "Size of signal must be the multiple of 2");

  if (N <= cutoff_size) {
	return dft(signal);
  }

  std::vector<double> even_signal;
  std::vector<double> odd_signal;

  even_signal.reserve(N / 2);
  odd_signal.reserve(N / 2);

  for (size_t i = 0; i < N; i++) {
	if (i % 2 == 0) {
	  even_signal.push_back(signal[i]);
	} else {
	  odd_signal.push_back(signal[i]);
	}
  }

  auto even = fft(even_signal);
  auto odd = fft(odd_signal);

  std::vector<std::complex<double>> factor_beg, factor_end;
  factor_beg.reserve(N / 2);
  factor_end.reserve(N / 2);

  for (size_t i = 0; i < N; i++) {
	auto val = (-2.0i * kPi * (double)i) / (double)N;
	if (i < N / 2) {
	  factor_beg.emplace_back(val);
	} else {
	  factor_end.emplace_back(val);
	}
  }

  auto result = concat<std::complex<double>>(
	  {add(even, mul(factor_beg, odd)), add(even, mul(factor_end, odd))});
  //  for (auto& x : result) { x = x * 1.0 / (double)N; }
  return result;
}
#endif// FOURIER_WITH_FFT
	
std::complex<double> dft2d(const std::vector<std::vector<double>> &signal,
                           int x, int y) {

  using namespace std::complex_literals;
  int N = f.size();
  int M = f.begin()->size();

  std::complex<double> sum;
  for (int k = 0; k < N; k++) {
    for (int l = 0; l < M; l++) {
      sum += f[k][l] * exp((-1i * 2.0 * PI) * (((double)k * (double)x) / N +
                                               ((double)l * (double)y) / M));
    }
  }
  return (1.0 / N * M) * sum;
}

std::vector<std::vector<double>>
dft2d(const std::vector<std::vector<double>> &signal) {
  using namespace std::complex_literals;
  std::vector<std::vector<double>> solution;
  for (int i = 0; i < f.size(); i++) {
    solution.emplace_back();
    for (int j = 0; j < f.begin()->size(); ++j)
      solution.back().emplace_back(0.0, 0.0);
  }

  int N = f.size();
  int M = f.begin()->size();

  for (int n = 0; n < N; n++) {
    for (int m = 0; m < M; m++) {
      std::complex<double> sum(0, 0);
      for (int k = 0; k < N; k++) {
        for (int l = 0; l < M; l++) {
          sum +=
              f[k][l] * exp((-1i * 2.0 * PI) * (((double)k * (double)n) / N +
                                                ((double)l * (double)m) / M));
        }
      }
      solution[n][m] = (1.0 / N * M) * sum;
    }
  }

  return solution;
}

}// namespace fourier
#endif

#endif// FOURIER_TRANSFORM_SRC_FOURIER_H_
