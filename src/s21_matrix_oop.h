#ifndef S21_MATRIX_OOP_H
#define S21_MATRIX_OOP_H

#include <initializer_list>
#include <iostream>

class S21Matrix {
 public:
  S21Matrix(const int rows = 1, const int cols = 1);

  S21Matrix(const S21Matrix& matrix);
  S21Matrix(S21Matrix&& matrix) noexcept;
  S21Matrix(const std::initializer_list<std::initializer_list<double>>& list);

  ~S21Matrix();

  int GetRows() const;
  int GetCols() const;
  void SetRows(const int rows);
  void SetCols(const int cols);

  bool EqMatrix(const S21Matrix& matrix) const;
  S21Matrix& SumMatrix(const S21Matrix& matrix);
  S21Matrix& SubMatrix(const S21Matrix& matrix);
  S21Matrix& MulNumber(const double num);
  S21Matrix& MulMatrix(const S21Matrix& matrix);
  S21Matrix Transpose() const;
  S21Matrix CalcComplements() const;
  double Determinant() const;
  S21Matrix InverseMatrix() const;

  S21Matrix& operator=(const S21Matrix& matrix);
  S21Matrix& operator=(S21Matrix&& matrix) noexcept;

  bool operator==(const S21Matrix& matrix) const;
  S21Matrix& operator+=(const S21Matrix& matrix);
  S21Matrix& operator-=(const S21Matrix& matrix);
  S21Matrix& operator*=(const S21Matrix& matrix);
  S21Matrix& operator*=(const double num);

  double& operator()(const int i, const int j);
  const double& operator()(const int i, const int j) const;

  friend S21Matrix operator+(const S21Matrix& a, const S21Matrix& b);
  friend S21Matrix operator-(const S21Matrix& a, const S21Matrix& b);
  friend S21Matrix operator*(const S21Matrix& a, const S21Matrix& b);
  friend S21Matrix operator*(const S21Matrix& matrix, const double num);
  friend S21Matrix operator*(const double num, const S21Matrix& matrix);

  // friend std::ostream& operator<<(std::ostream& out, const S21Matrix&
  // matrix);

 private:
  int rows_;
  int cols_;
  double** matrix_ = nullptr;

  void Create();
  void Copy(const S21Matrix& matrix);
  void Copy(const S21Matrix& matrix, const int rows, const int cols);
  void Move(S21Matrix&& matrix);
  void Clear();

  double GetMinor(const int row, const int col) const;

  void ThrowOutOfRange(const int i, const int j) const;
  void ThrowIsNotSquare() const;
  void ThrowDifferentDimensions(const S21Matrix& matrix) const;
};

#endif
