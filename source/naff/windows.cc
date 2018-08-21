#include "windows.h"

std::vector<std::complex<double>> cheb_window(const size_t N, const double a = 5.0) {
  std::vector<double> out;
  std::vector<std::complex<double>> complex_out;
  out.resize(N);
  auto cheb_poly = [](int n, double x) {
    double res;
    if (fabs(x) <= 1) res = cos(n*acos(x));
    else  res = cosh(n*acosh(x));
    return res;
  };
  double max=0;
  double tg = pow(10, a);
  double x0 = cosh((1.0/(N-1))*acosh(tg));
  double M = (N-1)/2;
  if(N%2==0) M = M + 0.5; 
  for(size_t nn=0; nn<(N/2+1); nn++){
      double n = nn-M;
      double sum = 0;
      for(size_t i=1; i<=M; i++){
          sum += cheb_poly(N-1,x0*cos(pi*i/N))*cos(2.0*n*pi*i/N);
      }
      out[nn] = tg + 2*sum;
      out[N-nn-1] = out[nn];
      if (out[nn]>max)
        max=out[nn];
  }
  for(size_t nn=0; nn<N; nn++) 
    out[nn] /= max; 
  for (const auto i:out) {
    complex_out.emplace_back(i,0.);
  } 
  return complex_out;
}

std::vector<std::complex<double>> hann_harm_window(const size_t N, const double n = 3.0) {
  std::vector<std::complex<double>> out;
  std::vector<std::complex<double>> out2;
  double T1 = 0;
  double T2 = N;
  double TM = (T2-T1)/2.0;
  double PIST = pi/TM;
  int factorial_1 = 1, factorial_2 = 1;
  for (size_t j=1;j<n+1;j++) {
    factorial_1 *= j;
  }
  for (size_t j=1;j<2*n+1;j++) {
    factorial_2 *= j;
  }
  double cn = pow(2, n) * pow(factorial_1, 2)/(factorial_2*1.0);
  for (size_t i=1; i<N+1; i++ ) {
    double T = (i-1) - TM;
    out.emplace_back(cn*pow((1.0 + cos(T*PIST)), n));
  }
  for (size_t i=0;i<N;i++) {
    const double multiplier = (1-cos(2*pi*i/(N-1.0)));
    out2.emplace_back(multiplier,0);
  }
  return out;
}

void WindowFunc::compute(const size_t N) {
  if (type == 'c') { 
    window = cheb_window(N, parameter);
  }
  else if (type == 'h') {
    window = hann_harm_window(N, parameter);
  }
  else if (type == 't') {
    window = taylor_window(N, parameter);
  }
  else if (type == 'n') {
    window = no_window(N, parameter);
  }
  else {
    throw std::runtime_error("window type is not defined");
  }
}

void WindowFunc::compute(const size_t N, const double param, const char tp) {
  if (type == 'c') { 
    window = cheb_window(N, param);
  }
  else if (type == 'h') {
    window = hann_harm_window(N, param);
  }
  else if (type == 't') {
    window = taylor_window(N, param);
  }
  else if (type == 'n') {
    window = no_window(N, param);
  }
  else {
    throw std::runtime_error("window type is not defined");
  }
   parameter = param;
   type = tp;
}

std::complex<double> WindowFunc::operator()(size_t t, size_t N) {
  if (window.size()!=N) {
    WindowFunc::compute(N);
  }
  return window[t];
}
 
std::complex<double> WindowFunc::operator()(size_t t, size_t N) const {
  if (window.size()!=N) {
    throw std::runtime_error("window and signal sizes do not match"); 
  }
  return window[t];
} 

std::vector<std::complex<double>> taylor_window(const size_t N, const double SLL) {
  std::vector<std::complex<double>> out;
  std::vector<double> ww(N);
  double A = acosh(pow(10, SLL)) / pi;
  double NBAR = round((2 * pow(A, 2) + 0.5)+0.5);

  //auto calculateFm = [] (int m, double sp2, double a, int nBar) {
  auto calculateFm = [] (size_t m, double sp2, double a, size_t nBar) {
    std::vector<std::vector<int>> n(nBar-1,std::vector<int>(1,0));
    std::vector<std::vector<int>> p(nBar-1,std::vector<int>(1,0));
    //int n[nBar - 1][1], p[nBar - 1][1];
    std::vector<std::vector<double>> Num(nBar-1,std::vector<double>(1,0));
    std::vector<std::vector<double>> Den(nBar-1,std::vector<double>(1,0));
    std::vector<std::vector<double>> Fm(nBar-1,std::vector<double>(1,0));
    double prodNum = 1, prodDen = 1, prodFm = 0.0;
    //double Num[nBar - 1][1], Den[nBar - 1][1], Fm[nBar - 1][1], prodNum = 1, prodDen = 1, prodFm = 0.0;
    for (size_t i = 0; i < nBar - 1; i++) {
      for (size_t j = 0; j < 1; j++) {
        n[i][j] = i + 1;
        p[i][j] = i + 1;
        Fm[i][j] = 0;
        Num[i][j] = 1 - (((m * m) / sp2)) / (a * a + (n[i][j] - 0.5) * ((n[i][j]) - 0.5));
        prodNum = prodNum * Num[i][j];
        Den[i][j] = 1 - (m * m / (double) (p[i][j] * p[i][j]));
        if ((i + 1) != m)  {
          prodDen = prodDen * Den[i][j];
        }
      }
    }
    prodFm = ((pow(-1, (m + 1)) * prodNum) / (2 * prodDen));
    return prodFm;
  };

  auto taylorWindow = [&calculateFm] (const size_t& N,double& NBAR, const double& SLL)  {
    double A, SP2, summation = 0.0;
    std::vector<double> Fm(NBAR);
    std::vector<size_t> k(N);
    std::vector<double> W(N);
    std::vector<double> sum(N);
    std::vector<double> Xi(N);
    A = acosh(pow(10, SLL)) / pi;
    SP2 = (NBAR * NBAR) / ((A * A + (NBAR - 0.5) * (NBAR - 0.5)));
    if (NBAR < (2 * pow(A, 2) + 0.5))  {
      std::cout<<"Warning: SSL value does not satisfy NBAR >= 2*A^2+0.5\n";
      std::cout<<" A = "<<A<<std::endl;
      std::cout<<" SLL = "<<SLL<<std::endl;
      std::cout<<" NBAR = "<<NBAR<<std::endl;
      std::cout<<" NBAR should be >= "<<2.0*pow(A,2)+0.5<<std::endl;
    }
    for (size_t i = 0; i < N; i++) {
      W[i] = 1.0;
      k[i] = i;
      Xi[i] = (k[i] - 0.5 * N + 0.5) / N;
    }
    for (size_t i = 0; i < NBAR - 1; i++) {
      Fm[i] = calculateFm((i + 1), SP2, A, NBAR);
    }
    for (size_t i = 0; i < N; i++) {
      summation = 0.0;
      for (size_t j = 1; j < NBAR; j++) {
        summation += Fm[j - 1] * cos(2 * pi * j * Xi[i]);
      }
      sum[i] = summation;
      W[i] = W[i] + 2 * sum[i];
    } 
    return W;
  };

  ww = taylorWindow(N, NBAR, SLL);
  for (size_t i = 0; i < N; i++){
    out.emplace_back(ww[i],0.);
  }
  return out;
}

std::vector<std::complex<double>> no_window(const size_t N, const double param) {
  std::vector<std::complex<double>> out;
  std::string message = "No window " +std::to_string(param);
  Print_opt::Write(Print_opt::Debug,message);
  for (size_t i = 0; i < N; i++){
    out.emplace_back(1.,0.);
  }
  return out;
}

