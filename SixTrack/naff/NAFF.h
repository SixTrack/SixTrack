#pragma once


#include <fftw3.h>

#include "signal.h"
#include "windows.h"
#include "spline_interpolation.h"
#include "minima.hpp"

typedef double Tfloat;
const std::complex<double> z (0,1);
typedef std::vector<double> double_vec;
typedef std::vector<std::complex<double>> complex_vec;

class NAFF {
  private:
    
  WindowFunc window;
  fftw_plan fftw_plan_;
  double_vec frequencies;
  double_vec amplitudes;
  size_t fft_size, f_counter;
  Signal signal, signal_no_upsampling;
  double max_index, fft_frequency = -1;
  std::vector<ComponentVector> norm_vectors;
  std::string merit_func, upsampling_type = "spline";
  bool f_found = true, flag_upsampling =false, flag_interpolation = false;
  double min_frequency =0;
  double max_frequency=1.0;
  double upsampling_factor = 10.0; 
  //////// Initialize signal and window
  double_vec remove_dc_component(double_vec &data) {
    double sum  = std::accumulate(data.begin(), data.end(), 0.0);
    double mean = sum / data.size();
    std::transform(data.begin(), data.end(), data.begin(), [mean](double x) { return x - mean; });
    return data;
  }
  void input(double_vec &init_data_x, double_vec &init_data_xp) {
    size_t new_size;
    remove_dc_component(init_data_x);
    remove_dc_component(init_data_xp);
    if (flag_interpolation == true) {
      new_size = multiple_of_six(init_data_x);
      Print_opt::Write(Print_opt::Debug,"----------> Hardy's interpolation used ...");
    }
    else {
      new_size = init_data_x.size();
    }
    if (flag_upsampling == true) {
      Print_opt::Write(Print_opt::Debug,"----------> Upsampling ON ...");
      for (size_t i = 0; i<new_size; i++) { 
        signal_no_upsampling.data.emplace_back(std::complex<double>(init_data_x[i], init_data_xp[i]));
      }
      if (upsampling_type == "spline") {
        for (double i = 1; i<new_size-1; i+=1.0/upsampling_factor) {
          signal.data.emplace_back(signal_no_upsampling[i]); 
        }
      }
      else if (upsampling_type == "linear") {
        for (double i = 1; i<new_size-1; i+=1.0/upsampling_factor) {
          signal.data.emplace_back(signal_no_upsampling(i)); 
        }
      }
      else 
        throw std::runtime_error("Interpolation method not defined! Options are linear or spline ");
      window.compute(signal.size());
    }
    else {
      for (size_t i=0; i<new_size; i++) { 
        signal.data.emplace_back(std::complex<double>(init_data_x[i], init_data_xp[i]));
      }
      window.compute(signal.size());
    }
  }
  //////// Fast Fourier Transform
  void FFTw () {
    std::vector<std::pair<double, double>> fftw_(signal.size());
    fftw_plan_ = fftw_plan_dft_1d(signal.size(), 
                   reinterpret_cast<fftw_complex*>(&signal.data[0]), 
                   reinterpret_cast<fftw_complex*>(fftw_.data()),
                   FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(fftw_plan_);
    fft_size = fftw_.size();
    max_fft_frequency(fftw_);
  }
  
  //////// First estimation of the peak frequency from FFT 
  void max_fft_frequency (std::vector<std::pair<double,double>> &fftw_) {
    double_vec amps;
    double_vec all_freqs;
    double max_amplitude = sqrt(fftw_[0].first*fftw_[0].first+fftw_[0].second*fftw_[0].second);
    max_index = 0.0;
    for (size_t i=0; i<fft_size; i++ ) {
      const double current_amplitude = sqrt(fftw_[i].first*fftw_[i].first+fftw_[i].second*fftw_[i].second);
      amps.push_back(current_amplitude);
      if (i<=fft_size/2)
        all_freqs.push_back(((i*1.0)/ (fft_size*1.0)));
      else
        all_freqs.push_back(1.0-(((i*1.0)/ (fft_size*1.0))));
      if (current_amplitude > max_amplitude && (abs(all_freqs.back())>min_frequency) && (abs(all_freqs.back()))<max_frequency) {
        max_amplitude = current_amplitude;
	max_index = i;
      }      
    }
   if (max_index<=fft_size/2)
      fft_frequency = ((max_index*1.0)/ (fft_size*1.0));
    else
      fft_frequency = (((max_index*1.0)/ (fft_size*1.0)-1.0));
  }
  
  //////// Maximization of <f(t),exp(i*2*pi*f*t)> function for refined frequency f
  ///////////////////// First method: Golden Section Search
  template <typename FuncMin>
  double cmp_min (const FuncMin& f, double min_value, double med_value,double max_value, double precision) {
    double phi = (1+sqrt(5))/2;
    double resphi=2-phi;
    if (std::abs(min_value-max_value)<precision) {
      return (min_value+max_value)/2;
    }
    double d = med_value+resphi*(max_value-med_value);
    if (f(d) < f(med_value))
          return cmp_min(f, med_value, d, max_value,precision);
    else
          return cmp_min(f, d, med_value, min_value, precision);
   }

  ///////////////////// Second method: BOOST Brent
  double maximize_fourier_integral () {
    auto y = [this](double f) { 
      Component c(f,signal.size());
      return (-abs(inner_product(signal, c, window, flag_interpolation))) ;};
    double step = 1.0/(signal.size());
    double min = fft_frequency - 1.0*step;
    double max = fft_frequency + 1.0*step;
    int bits = std::numeric_limits<double>::digits;
    std::pair<double, double> r = brent_minima(y, min, max, bits);
    Print_opt::Write(Print_opt::Debug,"Merit function: Maximize fourier integral");
    return r.first;   
  }
  
  ///////////////////// Subtraction of frequency components from signal 
  void subtract_frequency(Signal& signal, double& frequency) {
    Component v_i_pos(frequency, signal.size());
    Component v_i_neg(-frequency, signal.size());
    ComponentVector u_i(v_i_pos);
    if (f_counter!= 0) {
      for (size_t i=0; i<f_counter; i++) {
	u_i -= projection (v_i_pos, norm_vectors[i], window, flag_interpolation);
	u_i -= projection (v_i_neg, norm_vectors[i], window, flag_interpolation);
      } 
    }
    signal_projection (signal, u_i, window, flag_interpolation);
    signal -= u_i;
    amplitudes.push_back(abs(u_i.ampl)*signal.size());
    norm_vectors.push_back(u_i);
  }

  ///////////////////// Modified Gram Schmidt 
  /*void subtract_frequency(Signal& signal, double& frequency) {
    Component v_i(frequency, signal.size());
    ComponentVector u_i(v_i);
    std::vector<ComponentVector> u_k;
    u_k.emplace_back(u_i);
    if (f_counter!= 0) {
      for (size_t i=0; i<f_counter; i++) {
            u_i -= projection (u_k.back(), norm_vectors[i], window);
            u_k.emplace_back(u_i);
      }
    }
    signal_projection (signal, u_i, window);
    signal -= u_i;
    norm_vectors.push_back(u_i);
  }*/

  ////////////////////// Keep frequency which results in the minimum RMS in time domain
  double minimize_RMS_time () {
    auto y = [this](double current_frequency) {
      Signal signal_copy = signal;    
      Component v_i(current_frequency, signal_copy.size());
      subtract_frequency(signal_copy, current_frequency);
      double sum = 0;
      for (size_t i = 0; i <= signal_copy.size(); i++) {
        auto curr = abs(signal_copy[i]);
        sum += pow((curr),2);
       }
      return sqrt(sum);
    };   
    double step = 1.0/signal.size();
    double min = fft_frequency - step;
    double max = fft_frequency + step;
    int bits = std::numeric_limits<double>::digits;
    std::pair<double, double> r = brent_minima(y, min, max, bits);
    Print_opt::Write(Print_opt::Debug,"Merit function: Time domain energy");
    return r.first;
  }

  ////////////////////// Keep frequency which results in the minimum RMS in frequency domain
  double minimize_RMS_frequency () {
    int counter_now=0;
    auto y = [this, &counter_now](double current_frequency) {
      Signal signal_copy = signal;    
      Component v_i(current_frequency, signal_copy.size());
      subtract_frequency(signal_copy, current_frequency);
      double step = 1.0/signal.size();
      double min = fft_frequency - 1.0*step;
      double max = fft_frequency + 1.0*step;
      double stepp = (max-min)/100.0;

      auto fourier_integral = [&signal_copy, this](double f) { 
        Component c(f,signal_copy.size());
	return (-abs(inner_product(signal_copy, c,window, flag_interpolation))); };
      
      auto area= [this, &stepp, &max, &min, &fourier_integral, &current_frequency] () {
        double sum = 0;
        for (double i = min+stepp; i <= max; i+=stepp) {
	  auto curr = fourier_integral(i);
	  sum += pow((curr),2);
        }
        return sum;
      };
      return area();
    };
      
    double step = 1.0/signal.size();
    double min = fft_frequency - step;
    double max = fft_frequency + step;
    int bits = std::numeric_limits<double>::digits;
    std::pair<double, double> r = brent_minima(y, min, max, bits);
    Print_opt::Write(Print_opt::Debug,"Merit function: Fourier integral after subtraction");
    return r.first;
  }

  public:
  size_t fmax = 4;

  ~NAFF() { 
    fftw_destroy_plan(fftw_plan_);
  } 
  
  void set_window_parameter(const double p, const char tp) {
    window.parameter = p;
    window.type = tp;
  }

  double get_window_parameter() const {
    return window.parameter;
  }

  void set_merit_function(const std::string m) {
    merit_func = m;    
  }
  
  void set_upsampling(const bool flag, const double ups = 10.0, const std::string tp="spline") {
    flag_upsampling = flag;    
    upsampling_type = tp;
    upsampling_factor = ups;
  }
  
  void set_interpolation(const bool flag) {
    flag_interpolation = flag;    
  }

  void set_frequency_interval(const double& min_freq, const double& max_freq) {
    min_frequency = min_freq;
    max_frequency = max_freq;    
  }

  double_vec return_amplitudes() {
    return amplitudes;
  }
  
  //double_vec get_f1 (double_vec &init_data_x,double_vec &init_data_xp) {
  double get_f1 (double_vec &init_data_x,double_vec &init_data_xp) {
    Print_opt::SetLevel(2);
    if ( init_data_x.size()==0 && init_data_xp.size()==0 ) {
      fprintf(stderr, "CRITICAL ERROR in NAFF.h: x and xp input arrays are both empty");
      exit(EXIT_FAILURE);
    }
    if (frequencies.size() == 0) {
      input(init_data_x, init_data_xp);
    }
    FFTw();
    if (f_found == true) { 
      if (merit_func == "minimize_RMS_frequency") {
        frequencies.push_back(minimize_RMS_frequency());
      }
      else if (merit_func == "minimize_RMS_time") {
        frequencies.push_back(minimize_RMS_time());
      }
      ////////Default: maximize fourier integral
      else 
        frequencies.push_back(maximize_fourier_integral());
    }
    if (flag_upsampling == true) {
      for (auto& i:frequencies) {
        i*=upsampling_factor;
      }
    }
      return abs(frequencies.back());
  }
  
  double_vec get_f(double_vec &init_data_x, double_vec &init_data_xp) {
    f_counter = 0;
    Print_opt::SetLevel(2);
    while ((f_counter<fmax) && f_found == true) {
      std::string message = "Frequency: " + std::to_string(f_counter+1);	    
      Print_opt::Write(Print_opt::Debug, message);
      get_f1(init_data_x, init_data_xp);
      if (f_found == true) {
        subtract_frequency(signal, frequencies.back());
      }
    f_counter++;    
    }    
    std::string message = "Total number of frequencies found: "+std::to_string(f_counter); 
    Print_opt::Write(Print_opt::Debug, message);
    double_vec abs_frequencies;
    for (auto i:frequencies)
      abs_frequencies.push_back(std::abs(i));
    return abs_frequencies;
  }
};

