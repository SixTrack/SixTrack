#include "signal.h"
#include "windows.h"

template <typename T1, typename T2>
std::complex<double> inner_product(const T1& t1, const T2& t2, const WindowFunc & window_func, const bool& interpolation) {
  if ( t1.size() != t2.size()) throw std::runtime_error("Inner product sizes not equal!");
  if (interpolation == false) {
    std::complex<double> result;
    for (size_t i =0; i<t1.size(); i++) {
      result += t1[i]*std::conj(t2[i])*window_func(i,t1.size());
    }
    return result/(t1.size()*1.0);
  }
  else {
    std::complex<double> result2;
    std::vector<std::complex<double>> hardy_result;
    for (size_t i =0; i<t1.size(); i++) {
      hardy_result.emplace_back(t1[i]*std::conj(t2[i])*window_func(i,t1.size()));   
    }
    bool flag = false;
    size_t new_size = t1.size();
    while (flag==false) {
      if (new_size%6 == 1)
        flag = true;
      else
        new_size-=1;
    }
    size_t K = (new_size-1)/6; 
    result2 = 41.*hardy_result[0] + 216.*hardy_result[1]+27.* hardy_result[2] + 272.*hardy_result[3]+ 27.*hardy_result[4]+ 216.*hardy_result[5]+ 41.*hardy_result[new_size-1];
    for (size_t i=1;i<K;i++) {
      result2 += 82.0*hardy_result[6*i+1-1] + 216.0*hardy_result[6*i+2-1]+ 27.0*hardy_result[6*i+3-1]+272.0*hardy_result[6*i+4-1]+27.0*hardy_result[6*i+5-1]+216.0*hardy_result[6*i+6-1];
    }  
    double h = 1.0/(new_size*1.0);
    result2 *= h*6.0/840.0;
    return result2;  
  }
}

template std::complex<double>  inner_product<Component, ComponentVector>(const Component&, const ComponentVector&, const WindowFunc&, const bool&);
template std::complex<double>  inner_product<Signal, Component>(const Signal&, const Component&, const WindowFunc&, const bool&);
template std::complex<double>  inner_product<Signal, ComponentVector>(const Signal&, const ComponentVector&, const WindowFunc&, const bool&);

ComponentVector projection (const Component& v, const ComponentVector& u, const WindowFunc & window_func, const bool& interpolation){
  std::complex<double> num = inner_product(v, u, window_func, interpolation);
  std::complex<double> den = inner_product(u, u ,window_func, interpolation);
  ComponentVector s = v;
  s.ampl *= num/den;
  return s;
}

ComponentVector projection (const ComponentVector& v, const ComponentVector& u, const WindowFunc & window_func, const bool& interpolation){
  std::complex<double> num = inner_product(v, u, window_func, interpolation);
  std::complex<double> den = inner_product(u, u ,window_func, interpolation);
  ComponentVector s = v;
  s.ampl *= num/den;
  return s;
}

void signal_projection (const Signal& v, ComponentVector& u, const WindowFunc & window_func, const bool& interpolation){
  std::complex<double> num = inner_product(v, u, window_func, interpolation);
  std::complex<double> den = inner_product(u, u, window_func, interpolation);
  u.ampl = num/den;
}

double cmp_RMS (const std::vector<double>& data) {
  double res = abs(std::inner_product(data.begin(), data.end(), data.begin(),std::complex<double>(0,0))/(data.size()*1.0));
  return sqrt(res);  
}

std::ostream& operator <<(std::ostream& os, Component c) {
  os << "Frequency: " << c.freq << " Amplitude: " << c.ampl<< " ABS amplitude: "<< abs(c.ampl) << std::endl;
  return os;
}

std::ostream& operator <<(std::ostream& os, ComponentVector s) {
  for (size_t i=0;i<s.size();i++) {
    os << s[i] << std::endl;
  }
  return os;
}
std::ostream& operator <<(std::ostream& os, Signal s) {
  for (size_t i=0;i<s.size();i++) {
    os << s[i] << std::endl;
  }
  return os;
}

void write_file (const std::string& file_name, const std::vector<std::complex<double>>& var1 ,const double& min, const double& max, const double& step ) {
  std::ofstream myfile;
  myfile.open(file_name);
  for (double i=min; i<max; i+=step) {
    myfile<<std::setprecision(15)<<var1[i].real()<<" "<<std::setprecision(15)<<var1[i].imag()<<std::endl;
  }
}
void write_fft (const std::string& file_name, const std::vector<double>& var0,const std::vector<double>& var1 ,const double& min, const double& max, const double& step ) {
  std::ofstream myfile;
  myfile.open(file_name);
  for (double i=min; i<max; i+=step) {
    myfile<<std::setprecision(15)<<var0[i]<<" "<<var1[i]<<" "<<std::endl;
  }
}

void write_file_merit (const std::string& file_name, const std::function<double(double)>& y,const double& min, const double& max, const double& step ) {
  std::ofstream myfile;
  myfile.open(file_name);
  for (double i=min; i<max; i+=step) {
    myfile<<std::setprecision(14)<<i<<" "<<-y(i)<<std::endl;
  }
}

size_t multiple_of_six (std::vector<double>& data) {
  bool flag = false;
  size_t new_size = data.size();
  while (flag==false) {
    if (new_size%6 == 1)
      flag = true;
    else
      new_size-=1;
  }
  return new_size;
}

bool Print_opt::InitialisedM;
int Print_opt::levelM;

void Print_opt::Write(int level, std::string message){
  Initialised();
  if (level >= levelM) {
    std::cout<<message<<std::endl;
  }
}

void Print_opt::SetLevel(int level) {
  levelM = level;
  InitialisedM = true;
}

void Print_opt::Initialised() {
  if (!InitialisedM) {
    Init();
  }
}

void Print_opt::Init() {
  int default_level(Print_opt::All);
  SetLevel(default_level);
}

