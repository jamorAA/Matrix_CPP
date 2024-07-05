#include <gtest/gtest.h>

#include "s21_matrix_oop.h"

TEST(Test, DefaultConstructor) {
  S21Matrix matrix;
  ASSERT_EQ(matrix.GetRows(), 1);
  ASSERT_EQ(matrix.GetCols(), 1);
  ASSERT_EQ(matrix(0, 0), 0);
}

TEST(Test, Constructor1) {
  S21Matrix matrix(2, 4);
  ASSERT_EQ(matrix.GetRows(), 2);
  ASSERT_EQ(matrix.GetCols(), 4);
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 4; ++j) ASSERT_EQ(matrix(i, j), 0);
}

TEST(Test, Constructor2) {
  EXPECT_THROW(S21Matrix matrix(-2, 4), std::invalid_argument);
}

TEST(Test, CopyConstructor) {
  S21Matrix matrix = {{1, 2, 3.1}, {4.5, 5, 6.2}};
  S21Matrix copyMatrix = matrix;
  ASSERT_EQ(copyMatrix.GetRows(), 2);
  ASSERT_EQ(copyMatrix.GetCols(), 3);
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 3; ++j) ASSERT_EQ(matrix(i, j), copyMatrix(i, j));
}

TEST(Test, MoveConstructor) {
  S21Matrix matrix = {{1, 2, 3.1}, {4.5, 5, 6.2}};
  S21Matrix matrixTemp = matrix;
  S21Matrix copyMatrix = std::move(matrix);
  ASSERT_EQ(copyMatrix.GetRows(), 2);
  ASSERT_EQ(copyMatrix.GetCols(), 3);
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 3; ++j) ASSERT_EQ(matrixTemp(i, j), copyMatrix(i, j));
}

TEST(Test, EqMatrix1) {
  S21Matrix a = {{4.3, 3}, {11.3, -1.3}};
  S21Matrix b = {{4.3, 3}, {11.3, -1.3}};

  EXPECT_TRUE(a.EqMatrix(b));
}

TEST(Test, EqMatrix2) {
  S21Matrix a = {{4.3, 3}, {11.4, -1.3}};
  S21Matrix b = {{4.3, 3}, {11.3, -1.3}};

  EXPECT_FALSE(a.EqMatrix(b));
}

TEST(Test, EqMatrix3) {
  S21Matrix a = {{4.3, 3, 7}, {11.3, -1.3, -1}};
  S21Matrix b = {{4.3, 3}, {11.3, -1.3}};

  EXPECT_FALSE(a.EqMatrix(b));
}

TEST(Test, SumMatrix1) {
  S21Matrix a = {{4.3, 3, 7}, {11.3, -1.3, -1}};
  S21Matrix b = {{4.3, 3}, {11.3, -1.3}};

  EXPECT_THROW(a.SumMatrix(b), std::invalid_argument);
}

TEST(Test, SumMatrix2) {
  S21Matrix a = {{4, 3}, {1, -1}};
  S21Matrix b = {{9, 0}, {-5, 6}};

  S21Matrix res = {{13, 3}, {-4, 5}};

  EXPECT_TRUE(a.SumMatrix(b).EqMatrix(res));
}

TEST(Test, SumMatrix3) {
  S21Matrix a = {{4, 3}, {0, -1}};
  S21Matrix b = {{9, 0}, {-5, 6}};

  S21Matrix res = {{13, 3}, {-4, 5}};

  EXPECT_FALSE(a.SumMatrix(b).EqMatrix(res));
}

TEST(Test, SubMatrix1) {
  S21Matrix a = {{4.3, 3, 7}, {11.3, -1.3, -1}};
  S21Matrix b = {{4.3, 3}, {11.3, -1.3}};

  EXPECT_THROW(a.SumMatrix(b), std::invalid_argument);
}

TEST(Test, SubMatrix2) {
  S21Matrix a = {{4, 3}, {1, -1}};
  S21Matrix b = {{9, 0}, {-5, 6}};

  S21Matrix res = {{-5, 3}, {6, -7}};

  EXPECT_TRUE(a.SubMatrix(b).EqMatrix(res));
}

TEST(Test, SubMatrix3) {
  S21Matrix a = {{4, 3}, {0, -1}};
  S21Matrix b = {{9, 0}, {-5, 6}};

  S21Matrix res = {{-5, 3}, {6, -7}};

  EXPECT_FALSE(a.SubMatrix(b).EqMatrix(res));
}

TEST(Test, MulNumber) {
  S21Matrix a = {{4, 3}, {0, -1}};
  const double num = 2.5;

  S21Matrix res = {{10, 7.5}, {0, -2.5}};

  EXPECT_TRUE(a.MulNumber(num).EqMatrix(res));
}

TEST(Test, MulMatrix1) {
  S21Matrix a = {{1, 4}, {2, 5}, {3, 6}};
  S21Matrix b = {{1, -1, 1, 4}, {2, 3, 4, -5}, {1, 1, 1, 1}};

  EXPECT_THROW(a.MulMatrix(b), std::invalid_argument);
}

TEST(Test, MulMatrix2) {
  S21Matrix a = {{1, 4}, {2, 5}, {3, 6}};
  S21Matrix b = {{1, -1, 1}, {2, 3, 4}};

  S21Matrix res = {{9, 11, 17}, {12, 13, 22}, {15, 15, 27}};

  EXPECT_TRUE(a.MulMatrix(b).EqMatrix(res));
}

TEST(Test, Transponse) {
  S21Matrix a = {{1, 4}, {2, 5}, {3, 6}};

  S21Matrix res = {{1, 2, 3}, {4, 5, 6}};

  EXPECT_TRUE(res.EqMatrix(a.Transpose()));
}

TEST(Test, CalcComplements1) {
  S21Matrix a = {{1, 2, 3}, {0, 4, 2}, {5, 2, 1}};

  S21Matrix res = {{0, 10, -20}, {4, -14, 8}, {-8, -2, 4}};

  EXPECT_TRUE(res.EqMatrix(a.CalcComplements()));
}

TEST(Test, CalcComplements2) {
  S21Matrix a = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}, {1, 1, 1}};

  EXPECT_THROW(a.CalcComplements(), std::invalid_argument);
}

TEST(Test, Det1) {
  S21Matrix a = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}, {1, 1, 1}};

  EXPECT_THROW(a.Determinant(), std::invalid_argument);
}

TEST(Test, Det2) {
  S21Matrix a = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};

  EXPECT_EQ(a.Determinant(), 0);
}

TEST(Test, Det3) {
  S21Matrix a = {{1, 2, 3}, {4, 5, 6}, {7, 8, 10}};

  EXPECT_EQ(a.Determinant(), -3);
}

TEST(Test, Inverse1) {
  S21Matrix a = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};

  EXPECT_THROW(a.InverseMatrix(), std::invalid_argument);
}

TEST(Test, Inverse2) {
  S21Matrix a = {{2, 5, 7}, {6, 3, 4}, {5, -2, -3}};

  S21Matrix res = {{1, -1, 1}, {-38, 41, -34}, {27, -29, 24}};

  EXPECT_TRUE(res.EqMatrix(a.InverseMatrix()));
}

TEST(Test, Eq1) {
  S21Matrix a = {{4.3, 3}, {11.3, -1.3}};
  S21Matrix b = {{4.3, 3}, {11.3, -1.3}};

  EXPECT_TRUE(a == b);
}

TEST(Test, Eq2) {
  S21Matrix a = {{4.3, 3}, {11.5, -1.3}};
  S21Matrix b = {{4.3, 3}, {11.3, -1.3}};

  EXPECT_FALSE(a == b);
}

TEST(Test, Sum1) {
  S21Matrix a = {{4, 3}, {1, -1}};
  S21Matrix b = {{9, 0}, {-5, 6}};

  S21Matrix res = {{13, 3}, {-4, 5}};

  EXPECT_TRUE((a + b) == res);
}

TEST(Test, Sum2) {
  S21Matrix a = {{4, 3}, {9, -1}};
  S21Matrix b = {{9, 0}, {-5, 6}};

  S21Matrix res = {{13, 3}, {-4, 5}};

  EXPECT_FALSE((a + b) == res);
}

TEST(Test, Sum3) {
  S21Matrix a = {{4, 3}};
  S21Matrix b = {{9, 0}, {-5, 6}};

  EXPECT_THROW(a + b, std::invalid_argument);
}

TEST(Test, Sub1) {
  S21Matrix a = {{4, 3}, {1, -1}};
  S21Matrix b = {{9, 0}, {-5, 6}};

  S21Matrix res = {{-5, 3}, {6, -7}};

  EXPECT_TRUE((a - b) == res);
}

TEST(Test, Sub2) {
  S21Matrix a = {{4, 3}, {9, -1}};
  S21Matrix b = {{9, 0}, {-5, 6}};

  S21Matrix res = {{13, 3}, {-4, 5}};

  EXPECT_FALSE((a - b) == res);
}

TEST(Test, Sub3) {
  S21Matrix a = {{4, 3}};
  S21Matrix b = {{9, 0}, {-5, 6}};

  EXPECT_THROW(a - b, std::invalid_argument);
}

TEST(Test, MulNumberOp) {
  S21Matrix a = {{4, 3}, {0, -1}};
  const double num = 2.5;

  S21Matrix res = {{10, 7.5}, {0, -2.5}};

  EXPECT_TRUE(a * num == res);
}

TEST(Test, MulMatrixOp1) {
  S21Matrix a = {{1, 4}, {2, 5}, {3, 6}};
  S21Matrix b = {{1, -1, 1, 4}, {2, 3, 4, -5}, {1, 1, 1, 1}};

  EXPECT_THROW(a * b, std::invalid_argument);
}

TEST(Test, MulMatrixOp2) {
  S21Matrix a = {{1, 4}, {2, 5}, {3, 6}};
  S21Matrix b = {{1, -1, 1}, {2, 3, 4}};

  S21Matrix res = {{9, 11, 17}, {12, 13, 22}, {15, 15, 27}};

  EXPECT_TRUE(a * b == res);
}

TEST(Test, LeftSum1) {
  S21Matrix a = {{4, 3}, {1, -1}};
  S21Matrix b = {{9, 0}, {-5, 6}};

  S21Matrix res = {{13, 3}, {-4, 5}};

  EXPECT_TRUE((a += b) == res);
}

TEST(Test, LeftSum2) {
  S21Matrix a = {{4, 3}, {9, -1}};
  S21Matrix b = {{9, 0}, {-5, 6}};

  S21Matrix res = {{13, 3}, {-4, 5}};

  EXPECT_FALSE((a += b) == res);
}

TEST(Test, LeftSum3) {
  S21Matrix a = {{4, 3}};
  S21Matrix b = {{9, 0}, {-5, 6}};

  EXPECT_THROW(a += b, std::invalid_argument);
}

TEST(Test, LeftSub1) {
  S21Matrix a = {{4, 3}, {1, -1}};
  S21Matrix b = {{9, 0}, {-5, 6}};

  S21Matrix res = {{-5, 3}, {6, -7}};

  EXPECT_TRUE((a -= b) == res);
}

TEST(Test, LeftSub2) {
  S21Matrix a = {{4, 3}, {9, -1}};
  S21Matrix b = {{9, 0}, {-5, 6}};

  S21Matrix res = {{13, 3}, {-4, 5}};

  EXPECT_FALSE((a -= b) == res);
}

TEST(Test, LeftSub3) {
  S21Matrix a = {{4, 3}};
  S21Matrix b = {{9, 0}, {-5, 6}};

  EXPECT_THROW(a -= b, std::invalid_argument);
}

TEST(Test, LeftMulNumberOp) {
  S21Matrix a = {{4, 3}, {0, -1}};
  const double num = 2.5;

  S21Matrix res = {{10, 7.5}, {0, -2.5}};

  EXPECT_TRUE((a *= num) == res);
}

TEST(Test, LeftMulMatrixOp1) {
  S21Matrix a = {{1, 4}, {2, 5}, {3, 6}};
  S21Matrix b = {{1, -1, 1, 4}, {2, 3, 4, -5}, {1, 1, 1, 1}};

  EXPECT_THROW(a *= b, std::invalid_argument);
}

TEST(Test, LeftMulMatrixOp2) {
  S21Matrix a = {{1, 4}, {2, 5}, {3, 6}};
  S21Matrix b = {{1, -1, 1}, {2, 3, 4}};

  S21Matrix res = {{9, 11, 17}, {12, 13, 22}, {15, 15, 27}};

  EXPECT_TRUE((a *= b) == res);
}

TEST(Test, GetRows) {
  S21Matrix matrix(2, 2);
  EXPECT_TRUE(matrix.GetRows() == 2);
}

TEST(Test, GetCols) {
  S21Matrix matrix(2, 2);
  EXPECT_TRUE(matrix.GetRows() == 2);
}

TEST(Test, SetRows) {
  S21Matrix matrix(2, 2);
  matrix.SetRows(5);
  EXPECT_TRUE(matrix.GetRows() == 5);
}

TEST(Test, SetCols) {
  S21Matrix matrix(2, 2);
  matrix.SetCols(9);
  EXPECT_TRUE(matrix.GetCols() == 9);
}

TEST(Test, OutOfRange) {
  S21Matrix matrix = {{4, 5}, {3, 1}};
  EXPECT_THROW(matrix(1, 3), std::out_of_range);
}

int main(int argc, char *argv[]) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
