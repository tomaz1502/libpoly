#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <polyxx.h>

#include "doctest.h"

#include <random>

using namespace poly;

TEST_CASE("polynomial::properties") {
  Variable x("x");
  Variable y("y");
  Polynomial p1;
  Polynomial p2(Integer(5));
  Polynomial p3 = 3*x;
  Polynomial p4 = 3*x*y*y + 7*x;

  CHECK(is_zero(p1));
  CHECK_FALSE(is_zero(p2));
  CHECK_FALSE(is_zero(p3));
  CHECK_FALSE(is_zero(p4));

  CHECK(is_constant(p1));
  CHECK(is_constant(p2));
  CHECK_FALSE(is_constant(p3));
  CHECK_FALSE(is_constant(p4));

  CHECK_FALSE(is_linear(p1));
  CHECK_FALSE(is_linear(p2));
  CHECK(is_linear(p3));
  CHECK_FALSE(is_linear(p4));

  CHECK(is_lc_constant(p1));
  CHECK(is_lc_constant(p2));
  CHECK(is_lc_constant(p3));
  CHECK_FALSE(is_lc_constant(p4));

  CHECK(degree(p1) == 0);
  CHECK(degree(p2) == 0);
  CHECK(degree(p3) == 1);
  CHECK(degree(p4) == 2);

  CHECK(is_univariate(p3));
}

TEST_CASE("polynomial::resultant") {
  Variable y("y");
  Variable x("x");
  Polynomial p = 1 * pow(x, 6) + 2 * pow(x, 5) + 3 * y - 1;
  Polynomial q = 7 * pow(x, 5) + 5 * pow(x, 4);
  Polynomial r = resultant(p, q);
  CHECK(r == 28588707 * pow(y, 5) - 49925970 * pow(y, 4) + 34802730 * pow(y, 3) - 12107160 * pow(y, 2) + 2102235 * y - 145774);
}

TEST_CASE("polynomial::resultant_bezout") {
  Variable y("y");
  Variable x("x");

  SUBCASE("different degrees") {
    Polynomial p = 1 * pow(x, 6) + 2 * pow(x, 5) + 3 * y - 1;
    Polynomial q = 7 * pow(x, 5) + 5 * pow(x, 4);
    std::tuple<Polynomial, Polynomial, Polynomial> result = resultant_bezout(p, q);
    Polynomial res = std::get<0>(result);
    Polynomial u = std::get<1>(result);
    Polynomial v = std::get<2>(result);
    Polynomial expected_res = resultant(p, q);
    CHECK(res == expected_res);
    CHECK(u * p + v * q == res);
  }

  SUBCASE("equal degrees") {
    Polynomial p = pow(x, 2) + y * x + 1;
    Polynomial q = pow(x, 2) + x + y;
    std::tuple<Polynomial, Polynomial, Polynomial> result = resultant_bezout(p, q);
    Polynomial res = std::get<0>(result);
    Polynomial u = std::get<1>(result);
    Polynomial v = std::get<2>(result);
    Polynomial expected_res = resultant(p, q);
    CHECK(res == expected_res);
    CHECK(u * p + v * q == res);
  }

  SUBCASE("swapped order (deg(p) < deg(q))") {
    Polynomial p = 7 * pow(x, 5) + 5 * pow(x, 4);
    Polynomial q = 1 * pow(x, 6) + 2 * pow(x, 5) + 3 * y - 1;
    std::tuple<Polynomial, Polynomial, Polynomial> result = resultant_bezout(p, q);
    Polynomial res = std::get<0>(result);
    Polynomial u = std::get<1>(result);
    Polynomial v = std::get<2>(result);
    Polynomial expected_res = resultant(p, q);
    CHECK(res == expected_res);
    CHECK(u * p + v * q == res);
  }

  SUBCASE("linear polynomials") {
    Polynomial p = 2 * x + y;
    Polynomial q = 3 * x + 1;
    std::tuple<Polynomial, Polynomial, Polynomial> result = resultant_bezout(p, q);
    Polynomial res = std::get<0>(result);
    Polynomial u = std::get<1>(result);
    Polynomial v = std::get<2>(result);
    Polynomial expected_res = resultant(p, q);
    CHECK(res == expected_res);
    CHECK(u * p + v * q == res);
  }
}

TEST_CASE("polynomial::resultant_bezout_random") {
  // x must be declared last so it has the highest variable ID and becomes the
  // top variable in the polynomial ordering.
  Variable y("y");
  Variable z("z");
  Variable w("w");
  Variable t("t");
  Variable x("x");

  std::mt19937 rng(42);

  auto rand_int = [&](int lo, int hi) -> int {
    return std::uniform_int_distribution<int>(lo, hi)(rng);
  };

  auto rand_nonzero = [&]() -> int {
    int c = rand_int(1, 3);
    if (rand_int(0, 1)) c = -c;
    return c;
  };

  auto rand_poly = [&]() -> Polynomial {
    int deg = rand_int(1, 5);
    Polynomial result = rand_nonzero() * pow(x, deg);
    for (int d = deg - 1; d >= 0; --d) {
      int nterms = rand_int(0, 6);
      for (int j = 0; j < nterms; ++j) {
        Integer coeff(rand_nonzero());
        Polynomial term(coeff);
        if (rand_int(0, 1)) term = term * y;
        if (rand_int(0, 1)) term = term * z;
        if (rand_int(0, 1)) term = term * w;
        if (rand_int(0, 1)) term = term * t;
        if (d > 0) term = term * pow(x, d);
        result = result + term;
      }
    }
    return result;
  };

  for (int iter = 0; iter < 200; ++iter) {
    Polynomial p = rand_poly();
    Polynomial q = rand_poly();

    std::tuple<Polynomial, Polynomial, Polynomial> rb = resultant_bezout(p, q);
    Polynomial res = std::get<0>(rb);
    Polynomial u = std::get<1>(rb);
    Polynomial v = std::get<2>(rb);
    CHECK(res == resultant(p, q));
    CHECK(u * p + v * q == res);
  }
}

TEST_CASE("polynomial::discriminant") {
  Variable y("y");
  Variable x("x");
  Polynomial p = 1 * pow(x, 6) + 2 * pow(x, 5) + 3 * y - 1;
  Polynomial d = discriminant(p);
  CHECK(d == 11337408 * pow(y, 5) - 35095680 * pow(y, 4) + 34197120 * pow(y, 3) - 14999040 * pow(y, 2) + 3099840 * y - 246656);
}
TEST_CASE("polynomial::psc") {
  Variable y("y");
  Variable x("x");
  Polynomial p = 1 * pow(x, 6) + 2 * pow(x, 5) + 3 * y - 1;
  Polynomial q = 7 * pow(x, 5) + 5 * pow(x, 4);
  auto tmp = psc(p, q);

  CHECK(tmp.size() == 6);
  CHECK(tmp[0] == 28588707 * pow(y, 5) - 49925970 * pow(y, 4) + 34802730 * pow(y, 3) - 12107160 * pow(y, 2) + 2102235 * y - 145774);
  CHECK(tmp[1] == 416745 * y * y * y - 416745 * y * y + 138915 * y - 15435);
  CHECK(tmp[2] == Integer(0));
  CHECK(tmp[3] == Integer(0));
  CHECK(tmp[4] == Integer(-45));
  CHECK(tmp[5] == Integer(7));
}

TEST_CASE("polynomial::operator<<") {
  Variable y("y");
  Variable x("x");
  Polynomial p = 1 * pow(x, 6) + 2 * pow(x, 5) + 3 * y - 1;
  std::stringstream out;
  out << p;
  CHECK(out.str() == "1*x^6 + 2*x^5 + (3*y - 1)");
}

TEST_CASE("polynomial::constructor") {
  Variable y("y");
  Variable x("x");
  Polynomial p = 1 * pow(x, 6) + 2 * pow(x, 5) + 3 * y - 1;
  lp_polynomial_t* ptr = p.get_internal();

  Polynomial p_cp(p);
  Polynomial p_mv(std::move(p));
  CHECK(ptr == p_mv.get_internal());
  CHECK(ptr != p_cp.get_internal());
  CHECK(p_cp == p_mv);

  Polynomial p_mv_a, p_cp_a;
  p_cp_a = p_cp;
  p_mv_a = std::move(p_mv);
  CHECK(p_cp == p_cp_a);
  CHECK(ptr == p_mv_a.get_internal());
}

TEST_CASE("polynomial::evaluate") {
  Variable x("x");
  Variable y("y");
  Polynomial p1 = 1 * pow(x, 6) + 2 * pow(x, 5) + 3 * y - 1;
  Polynomial p2 = 1 * y * pow(x + 1, 4);

  // full assignment
  Assignment a;
  a.set(x, -1);
  a.set(y, 2);
  CHECK(a.has(x));
  CHECK(a.has(y));

  CHECK(evaluate(p1, a) == Value(4));
  CHECK(evaluate_constraint(p1, a, SignCondition::GE) == true);
  CHECK(evaluate_constraint(p1, a, SignCondition::LT) == false);
  CHECK(evaluate_constraint(p1, a, SignCondition::EQ) == false);
  CHECK(evaluate_constraint_subs(p1, a, SignCondition::GE) == 1);
  CHECK(evaluate_constraint_subs(p1, a, SignCondition::LT) == 0);
  CHECK(evaluate_constraint_subs(p1, a, SignCondition::EQ) == 0);

  CHECK(evaluate(p2, a) == Value(0));
  CHECK(evaluate_constraint(p2, a, SignCondition::GE) == true);
  CHECK(evaluate_constraint(p2, a, SignCondition::LT) == false);
  CHECK(evaluate_constraint(p2, a, SignCondition::EQ) == true);
  CHECK(evaluate_constraint_subs(p2, a, SignCondition::GE) == 1);
  CHECK(evaluate_constraint_subs(p2, a, SignCondition::LT) == 0);
  CHECK(evaluate_constraint_subs(p2, a, SignCondition::EQ) == 1);

  // partial assignment
  a.unset(y);
  CHECK(a.has(x));
  CHECK(!a.has(y));

  CHECK(evaluate_constraint_subs(p1, a, SignCondition::GE) == -1);
  CHECK(evaluate_constraint_subs(p1, a, SignCondition::LT) == -1);
  CHECK(evaluate_constraint_subs(p1, a, SignCondition::EQ) == -1);

  CHECK(evaluate_constraint(p2, a, SignCondition::GE) == true);
  CHECK(evaluate_constraint(p2, a, SignCondition::LT) == false);
  CHECK(evaluate_constraint(p2, a, SignCondition::EQ) == true);
  CHECK(evaluate_constraint_subs(p2, a, SignCondition::GE) == 1);
  CHECK(evaluate_constraint_subs(p2, a, SignCondition::LT) == 0);
  CHECK(evaluate_constraint_subs(p2, a, SignCondition::EQ) == 1);
}
