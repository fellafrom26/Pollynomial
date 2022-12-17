#include "polynomial.h"

bool is_prime(int num) {
	for (int i = 2; i * i <= num; ++i) {
		if (num % i == 0) {
			return false;
		}
	}

	return true;
}
int find_k(int x, int y) {
  int answer = 1;
  while ((answer * y) % p != x) {
    ++answer;
  }

  return answer;
}
std::vector<int64_t> get_divisors(int64_t num) {
  std::vector<int64_t> answer;
  if (num == 2 || num == 3) {
    answer.push_back(num);
    return answer;
  }

  int divisor = 2;
  while (divisor * divisor <= num) {
    if (num % divisor == 0) {
      answer.push_back(divisor);
      do {
        num /= divisor;
      } while (num % divisor == 0);
    } else if (divisor == 2) {
      divisor = 3;
    } else {
      divisor += 2;
    }
  }

  if (num != 1) {
    answer.push_back(num);
  }

  return answer;
}

Polynomial::Polynomial() : degree(0), body(1, 0) {
}
Polynomial::Polynomial(const int64_t &m) : degree(m) {
  body.resize(degree + 1, 0);
}
Polynomial::Polynomial(const int64_t &deg, const std::vector<int> &v) {
  if (deg == 0 && v[0] == 0) {
    degree = 0;
    body = std::vector<int>(1, 0);
    return;
  }
  int64_t i = deg;
  while (v[i] == 0 && i >= 0) {
    --i;
  }
  if (i == -1) {
    degree = 0;
    body = std::vector<int>(1, 0);
    return;
  }
  degree = i;
  body.resize(degree + 1);
  for (; i >= 0; --i) {
    body[i] = v[i];
  }
}
void Polynomial::normalize() {
  if (degree == 0) {
    body[0] %= p;
    return;
  }

  if (body[degree] != 0) {
    for (auto &it: body) {
      it %= p;
    }
    return;
  }

  int64_t i = degree;
  while (body[i] == 0 && i >= 0) {
    --i;
  }
  if (i == -1) {
    degree = 0;
    body = std::vector<int>(1, 0);
    return;
  }
  degree = i;

  std::vector<int> new_body;
  new_body.resize(degree + 1);
  while (i >= 0) {
    new_body[i] = body[i] % p;
    --i;
  }

  body = new_body;
}
int64_t Polynomial::get_degree() const {
  return degree;
}
int &Polynomial::operator[](int64_t index) {
  return body[index];
}
const int &Polynomial::operator[](int64_t index) const {
  return body[index];
}

Polynomial operator+(const Polynomial &a, const Polynomial &b) {
  if (b.get_degree() > a.get_degree()) {
    Polynomial answer = b;
    for (int64_t i = 0; i <= a.get_degree(); ++i) {
      answer[i] += a[i];
      if (answer[i] >= p) {
        answer[i] -= p;
      }
    }
    answer.normalize();
    return answer;
  }

  Polynomial answer = a;
  for (int64_t i = 0; i <= b.get_degree(); ++i) {
    answer[i] += b[i];
    if (answer[i] >= p) {
      answer[i] -= p;
    }
  }

  answer.normalize();
  return answer;
}
Polynomial operator-(const Polynomial &a, const Polynomial &b) {
  Polynomial answer = a;
  for (int64_t i = 0; i <= b.get_degree(); ++i) {
    answer[i] -= b[i];
    if (answer[i] < 0) {
      answer[i] += p;
    }
  }

  answer.normalize();
  return answer;
}
Polynomial operator*(const Polynomial &a, const Polynomial &b) {
  Polynomial answer(a.get_degree() + b.get_degree());

  for (int i = 0; i <= a.get_degree(); ++i) {
    for (int j = 0; j <= b.get_degree(); ++j) {
      answer[i + j] += (a[i] * b[j]) % p;
    }
  }

  answer.normalize();
  return answer;
}
Polynomial operator%(const Polynomial &a, const Polynomial &b) {
  if (b.get_degree() == 0) {
    return Polynomial{};
  }
  Polynomial x;
  Polynomial A(a);

  while (A.get_degree() >= b.get_degree()) {
    x = Polynomial{A.get_degree() - b.get_degree()};
    x[A.get_degree() - b.get_degree()] = find_k(A[A.get_degree()], b[b.get_degree()]);
    A -= x * b;
  }

  return A;
}
Polynomial &Polynomial::operator+=(const Polynomial &right) {
  *this = *this + right;
  return *this;
}
Polynomial &Polynomial::operator-=(const Polynomial &right) {
  *this = *this - right;
  return *this;
}
Polynomial &Polynomial::operator*=(const Polynomial &right) {
  *this = *this * right;
  return *this;
}
Polynomial &Polynomial::operator%=(const Polynomial &right) {
  *this = *this % right;
  return *this;
}

void Polynomial::pow_p() {
  degree *= p;
  std::vector<int> new_body(degree + 1, 0);

  for (int i = 0; i < body.size(); ++i) {
    new_body[i * p] = body[i];
  }

  body = new_body;
}
Polynomial gcd(const Polynomial &a, const Polynomial &b) {
  Polynomial temp, x(a), y(b);
  while (!(y.get_degree() == 0 && y[0] == 0)) {
    temp = y;
    y = x % y;
    x = temp;
  }
  return x;
}
bool Polynomial::is_irreducible() const {
  auto primes = get_divisors(degree);
  std::vector<int64_t> powers(primes.size());
  for (int64_t i = static_cast<int64_t>(primes.size()) - 1, j = 0; i >= 0; --i, ++j) {
    powers[j] = degree / primes[i];
  }
  int count = 0;

  Polynomial add(1);
  add[1] = p - 1;
  Polynomial g(1);
  g[1] = 1;
  g %= *this;

  for (int i = 1; i <= degree; ++i) {
    g.pow_p();
    g %= *this;

    if (i == powers[count]) {
      ++count;

      if (gcd((g + add), *this).get_degree() != 0) {
        return false;
      }
    }
  }

  g += add;
  g %= *this;
  if (g.get_degree() == 0 && g[0] == 0) {
    return true;
  }
  return false;
}
void Polynomial::print() const {
  int64_t i = degree;
  if (degree == 0) {
    std::cout << body[0] << '\n';
    return;
  }
  if (body[i] != 1) {
    std::cout << body[i];
  }
  std::cout << "x";
  if (i > 1) {
    std::cout << '^' << i;
  }
  --i;

  for (; i > 1; --i) {
    if (body[i] == 0) {
      continue;
    }
    if (body[i] > 0) {
      if (body[i] == 1) {
        std::cout << " + ";
      } else {
        std::cout << " + " << body[i];
      }
    }
    std::cout << "x^" << i;
  }
  if (body.size() > 2) {
    if (body[1] > 0) {
      std::cout << " + ";
      if (body[1] != 1) {
        std::cout << body[1];
      }
      std::cout << 'x';
    }
  }

  if (body[0] > 0) {
    std::cout << " + " << body[0];
  }
  std::cout << '\n';
}

int main() {
  srand(time(nullptr));

  clock_t begin, end;
  int ch;
  while (true) {
    std::cout << "Insert p: ";
    std::cin >> p;
	if (!is_prime(p)) {
		std::cout << "number is not prime!\n";
		continue;
	}
    std::cout << "Insert n: ";
    std::cin >> n;

    std::vector<int> temp(n + 1, 1);
    begin = clock() / CLOCKS_PER_SEC;
    for (int count = 1; true; ++count) {
      for (int i = n - 1; i >= 0; --i) {
        temp[i] = rand() % p;
      }
      Polynomial polynomial(n, temp);

      if (polynomial.get_degree() <= 1) {
        end = clock() / CLOCKS_PER_SEC;
        std::cout << "DONE! " << count << " attempts were made within ~" << (end - begin) << " s\n";
        polynomial.print();
        break;
      }
      if (polynomial.is_irreducible()) {
        end = clock() / CLOCKS_PER_SEC;
        std::cout << "DONE! " << count << " attempts were made within ~" << (end - begin) << " s\n";
        polynomial.print();
        break;
      }
    }

	std::cout << "\npress Enter to continue or type \"done\" to exit\n";
    getchar();
    ch = getchar();
    if (ch != '\n') {
      break;
    }
  }

  std::cout << "OK\n";
}
