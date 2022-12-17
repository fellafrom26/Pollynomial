#ifndef IRREDUCIBILITY_TEST__POLYNOMIAL_H_
#define IRREDUCIBILITY_TEST__POLYNOMIAL_H_

#include <iostream>
#include <vector>
#include <ctime>
#include <cstdint>

int p;
int n;

class Polynomial {
 private:
  std::vector<int> body;
  int64_t degree;

 public:
  Polynomial();
  explicit Polynomial(const int64_t &m);
  Polynomial(const int64_t &m, const std::vector<int> &v);
  void normalize();
  void print() const;
  void pow_p();
  int64_t get_degree() const;
  int &operator[](int64_t index);
  const int &operator[](int64_t index) const;
  bool is_irreducible() const;

  Polynomial &operator+=(const Polynomial &right);
  Polynomial &operator-=(const Polynomial &right);
  Polynomial &operator*=(const Polynomial &right);
  Polynomial &operator%=(const Polynomial &right);
};

#endif //IRREDUCIBILITY_TEST__POLYNOMIAL_H_