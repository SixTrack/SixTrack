#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <numeric>
#include <cmath>
#include <algorithm>
#include <functional>
#include <iomanip>

#include "consts.h"
#include "spline_interpolation.h"

class Component {
  public:
    double freq;
    size_t signal_size; 
    std::complex<double> ampl;
    
    Component (): ampl(1.,0.) {}    
    Component(double freq, size_t signal_size): freq(freq), signal_size(signal_size), ampl(1.,0.) {}
    ~Component () {}
    
    std::complex<double> operator[] (size_t t) const {
      return (exp(2*pi*freq*t*imag));    
    }
    
    size_t size() const {
      return signal_size;
    }
};

class ComponentVector { 
  public:
    typedef std::vector<std::complex<double>> data_t;
    std::complex<double> ampl;
    std::vector<Component> cs;

    ComponentVector():ampl(1,0) {}
    ComponentVector(const Component &other): ampl(1.,0.) { cs.emplace_back(other);}
    ~ComponentVector() {}

    size_t size() const{
      return cs[0].size();
    }

    std::complex<double> operator[] (size_t t) const {
      std::complex<double> res;
      for (const auto & c: cs) {
        res += ampl*c.ampl*c[t];
      }
      return res;
    }


    ComponentVector operator-=(const ComponentVector& other) {
      if (size() != other.size() ) throw std::runtime_error("Inner product sizes not equal!");
      for (auto i: other.cs) {
        i.ampl*=-other.ampl;
        cs.emplace_back(i);
      }
      return (*this);
    }
};

class Signal {
  public:
    ComponentVector::data_t data; 
    Signal() {}
    
    Signal(const std::vector<std::complex<double>> & v){
      for (const auto i:v) data.emplace_back(i);
    }
    Signal(const Signal &) = default;
    Signal& operator=(const Signal &) = default;

    ~Signal() {}
    
    std::complex<double> operator[] (size_t t) const {
      return data[t];
    }    

    std::complex<double> operator[] (double t) const {
      return spline(t, data);
     }

    std::complex<double> operator() (double t) const {
      size_t a = (int)(t);
      size_t b = (int)(t)+1;
      return data[a]+(data[b]-data[a])*(t-a)/(1.0*(b-a));
     }
           
    size_t size() const{
      return data.size();
    }

    Signal operator-=(const ComponentVector& other) {
      if (size() != other.size() ) throw std::runtime_error("Inner product sizes not equal!");
      for (size_t i = 0; i<size(); ++i) {
        data[i] -= other[i];
        data[i] -= std::conj(other[i]);
      }
      return *this;
    }    
};

class WindowFunc;

template <typename T1, typename T2>
std::complex<double> inner_product(const T1&, const T2&, const WindowFunc &, const bool&);

ComponentVector projection (const Component&, const ComponentVector &, const WindowFunc&, const bool&);
ComponentVector projection (const ComponentVector&, const ComponentVector&, const WindowFunc&, const bool&);

void signal_projection (const Signal&, ComponentVector&, const WindowFunc&, const bool&);

double cmp_RMS (const std::vector<double>&);

std::ostream& operator <<(std::ostream&, Signal);

std::ostream& operator <<(std::ostream& , ComponentVector);

std::ostream& operator <<(std::ostream&, Component);

void write_file (const std::string&, const std::vector<std::complex<double>>&, const double&, const double&, const double&) ;
void write_fft (const std::string&,const std::vector<double>&, const std::vector<double>&, const double&, const double&, const double&) ;

void write_file_merit (const std::string&, const std::function<double(double)>&, const double& , const double&, const double&);
  
size_t multiple_of_six(std::vector<double>&);

class Print_opt {
  public:
    enum {All = 0, Debug, Info};
    static void Write(int level, std::string message);
    static void SetLevel(int level);
  protected:
    static void Initialised();
    static void Init();
  private:
    Print_opt();
    static bool InitialisedM;
    static int levelM;
};


