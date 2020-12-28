#ifndef Binning_h
#define Binning_h

#include "TArrayD.h"
#include "TMath.h"  // Sort

//==============================================================================
// (1) Define the variable binning schemes for analysis variables
// (2) Some utility binning functions
//==============================================================================
TArrayD MakeUniformBinArray(int nbins, double min, double max);
void SortArray(TArrayD& arr);
TArrayD GetSortedArray(const TArrayD& arr);
TArrayD GetTArrayFromVec(const std::vector<double>& vec);

namespace CCPi {
TArrayD GetBinning(const std::string var_name) {
  std::vector<double> bins_vec;
  if (var_name == "enu") {
    bins_vec = {0., 1.e3, 3.e3, 4.e3, 6.5e3, 9.5e3, 14.e3, 30.e3};
  } else if (var_name == "pmu") {
    bins_vec = {0., 1.e3, 2.e3, 3.e3, 4.e3, 5.5e3, 7.5e3, 10.e3, 13.e3, 30.e3};
  } else if (var_name == "q2") {
    bins_vec = {0.e5, 1.e5, 2.e5, 4.e5, 7.e5, 10.e5, 13.e5, 20.e5};
  } else if (var_name == "thetamu_deg") {
    bins_vec = {0., 1., 2.,  3.,  4.,  5.,  6.,  7.,
                8., 9., 10., 11., 12., 14., 16., 20.};
  } else if (var_name == "thetapi_deg") {
    bins_vec = {0., 15., 30., 45, 60, 76., 108., 122., 136., 150., 165., 180.};
  } else if (var_name == "tpi" || var_name == "tpi_mbr") {
    bins_vec = {0., 35., 68., 100., 133., 166., 200., 350.};
  } else if (var_name == "wexp") {
    bins_vec = {10.e2, 11.e2, 12.e2, 13.e2, 14.e2, 15.e2};
  } else if (var_name == "ptmu") {
    bins_vec = {0.,   1.e2, 2.e2,  3.e2,   4.e2,  5.e2,
                6.e2, 8.e2, 10.e2, 12.5e2, 15.e2, 25.e2};
  } else if (var_name == "pzmu") {
    bins_vec = {0.,    1.e3, 2.e3,  3.e3,  4.e3, 5.e3,
                6.0e3, 8.e3, 10.e3, 15.e3, 20.e3};
  }
  // prepare an array from the bin vector
  TArrayD bins_array(GetTArrayFromVec(bins_vec));
  SortArray(bins_array);
  return bins_array;
}
}  // namespace CCPi

TArrayD MakeUniformBinArray(int nbins, double min, double max) {
  double step_size = (max - min) / nbins;
  double arr[nbins + 1];  // +1 because binning arrays include top edge.
  for (int i = 0; i <= nbins; ++i) arr[i] = min + i * step_size;
  const int size = sizeof(arr) / sizeof(*arr);
  return TArrayD(size, arr);
}

void SortArray(TArrayD& arr) {
  // store in index the correct order of arr
  int size = arr.GetSize();
  int* index = new int[size];
  TMath::Sort(size, arr.GetArray(), index,
              kFALSE);  // default ordering is decreasing. y tho.

  // then make a new array from that correct ordering
  TArrayD dummy(size);
  for (int i = 0; i < size; i++) dummy[i] = arr[index[i]];

  // finally set arr
  arr = dummy;
}

TArrayD GetSortedArray(const TArrayD& arr) {
  // store in index the correct order of arr
  int size = arr.GetSize();
  int* index = new int[size];
  TMath::Sort(size, arr.GetArray(), index,
              kFALSE);  // default ordering is decreasing. y tho.

  // then make a new array from that correct ordering
  double sorted_array[size];
  for (int i = 0; i < size; i++) sorted_array[i] = arr[index[i]];

  // finally return sorted
  return TArrayD(size, sorted_array);
}

TArrayD GetTArrayFromVec(const std::vector<double>& vec) {
  double array[vec.size()];
  std::copy(vec.begin(), vec.end(), array);
  const int size = sizeof(array) / sizeof(*array);
  return TArrayD(size, array);
}

#endif  // Binning_h