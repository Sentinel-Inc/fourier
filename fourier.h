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

///
/// \param signal
/// \return
double fft(const std::vector<std::complex<double>>& signal);

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

double fft(const std::vector<std::complex<double>>& signal) {

  assert(false && "Not implemented");
#if 0
   size_t N = signal.size();
   size_t cutoff_size = 32;

   // TODO: That check is not complete
   assert(N % 2 == 0 && "Size of signal must be the power of 2");

   if (N < cutoff_size){
	return dft(signal);
   } else {

	auto signal_even = fft();
	auto signal_odd = fft();
	auto factor = nullptr;
	return concat()
   }
#endif
  return -1;
}

}// namespace fourier
#endif

#endif// FOURIER_TRANSFORM_SRC_FOURIER_H_
