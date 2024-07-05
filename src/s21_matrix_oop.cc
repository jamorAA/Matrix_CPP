#include "s21_matrix_oop.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

void S21Matrix::Create() {
  matrix_ = new double*[rows_];
  for (int i = 0; i < rows_; ++i) matrix_[i] = new double[cols_]{0.0};
}

void S21Matrix::Copy(const S21Matrix& matrix) {
  for (int i = 0; i < matrix.rows_; ++i) {
    for (int j = 0; j < matrix.cols_; ++j) (*this)(i, j) = matrix(i, j);
  }
}

void S21Matrix::Copy(const S21Matrix& matrix, const int rows, const int cols) {
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) (*this)(i, j) = matrix(i, j);
  }
}

void S21Matrix::Move(S21Matrix&& matrix) {
  rows_ = matrix.rows_;
  matrix.rows_ = 0;

  cols_ = matrix.cols_;
  matrix.cols_ = 0;

  matrix_ = matrix.matrix_;
  matrix.matrix_ = nullptr;
}

void S21Matrix::Clear() {
  for (int i = 0; i < rows_; ++i) delete[] matrix_[i];
  delete[] matrix_;
}

// объединил два конструктора в один
S21Matrix::S21Matrix(const int rows, const int cols)
    : rows_(rows), cols_(cols) {
  if (rows_ < 1 || cols_ < 1)
    throw std::invalid_argument("negative or zero values for matrix dimension");

  Create();
}

S21Matrix::S21Matrix(const S21Matrix& matrix)
    : S21Matrix(matrix.rows_, matrix.cols_) {
  Copy(matrix);
}

S21Matrix::S21Matrix(S21Matrix&& matrix) noexcept { Move(std::move(matrix)); }

// добавил чисто ради удобства, знаю, что можно сломать
S21Matrix::S21Matrix(
    const std::initializer_list<std::initializer_list<double>>& list)
    : rows_(static_cast<int>(list.size())) {
  for (const auto& i : list) {
    cols_ = static_cast<int>(i.size());
    break;
  }

  Create();

  int i = 0, j = 0;
  for (const auto& k : list) {
    for (const auto& elem : k) matrix_[i][j++] = elem;
    j = 0;
    ++i;
  }
}

S21Matrix::~S21Matrix() { Clear(); }

int S21Matrix::GetRows() const { return rows_; }

int S21Matrix::GetCols() const { return cols_; }

void S21Matrix::SetRows(const int rows) {
  if (rows < 1)
    throw std::invalid_argument("negative or zero values for matrix dimension");

  S21Matrix temp = (*this);
  Clear();
  rows_ = rows;
  Create();
  Copy(temp, std::min(rows_, temp.rows_), std::min(cols_, temp.cols_));
}

void S21Matrix::SetCols(const int cols) {
  if (cols < 1)
    throw std::invalid_argument("negative or zero values for matrix dimension");

  S21Matrix temp = (*this);
  Clear();
  cols_ = cols;
  Create();
  Copy(temp, std::min(rows_, temp.rows_), std::min(cols_, temp.cols_));
}

bool S21Matrix::EqMatrix(const S21Matrix& matrix) const {
  if ((rows_ != matrix.rows_) || (cols_ != matrix.cols_)) return false;

  const double eps = 1e-7;
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      if (fabs((*this)(i, j) - matrix(i, j)) >= eps) return false;
    }
  }

  return true;
}

S21Matrix& S21Matrix::SumMatrix(const S21Matrix& matrix) {
  ThrowDifferentDimensions(matrix);

  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) (*this)(i, j) += matrix(i, j);
  }

  return *this;
}

S21Matrix& S21Matrix::SubMatrix(const S21Matrix& matrix) {
  ThrowDifferentDimensions(matrix);

  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) (*this)(i, j) -= matrix(i, j);
  }

  return *this;
}

S21Matrix& S21Matrix::MulNumber(const double num) {
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) (*this)(i, j) *= num;
  }

  return *this;
}

S21Matrix& S21Matrix::MulMatrix(const S21Matrix& matrix) {
  if (cols_ != matrix.rows_)
    throw std::invalid_argument(
        "the number of columns of the first matrix is not equal to the number "
        "of rows of the second matrix");

  S21Matrix temp = (*this);

  Clear();
  rows_ = temp.GetRows();
  cols_ = matrix.GetCols();
  Create();

  for (int i = 0; i < temp.rows_; ++i) {
    for (int j = 0; j < matrix.cols_; ++j) {
      double sum = 0;
      for (int k = 0; k < temp.cols_; ++k) sum += temp(i, k) * matrix(k, j);
      (*this)(i, j) = sum;
    }
  }

  return *this;
}

S21Matrix S21Matrix::Transpose() const {
  S21Matrix transMatrix(cols_, rows_);

  for (int i = 0; i < cols_; ++i) {
    for (int j = 0; j < rows_; ++j) transMatrix(i, j) = (*this)(j, i);
  }

  return transMatrix;
}

double S21Matrix::GetMinor(const int row, const int col) const {
  if (rows_ == 1) return 1.0;

  S21Matrix res(rows_ - 1, cols_ - 1);

  int iRes = 0, jRes = 0;

  for (int i = 0; i < rows_; ++i) {
    if (i == row) continue;

    for (int j = 0; j < cols_; ++j) {
      if (j != col) res(iRes, jRes++) = (*this)(i, j);
    }

    ++iRes;
    jRes = 0;
  }

  return res.Determinant();
}

S21Matrix S21Matrix::CalcComplements() const {
  ThrowIsNotSquare();

  S21Matrix res(rows_, cols_);

  for (int i = 0; i < res.rows_; ++i) {
    for (int j = 0; j < res.cols_; ++j)
      res(i, j) = GetMinor(i, j) * pow(-1, i + j);
  }

  return res;
}

double S21Matrix::Determinant() const {
  ThrowIsNotSquare();

  S21Matrix a = (*this);

  int index;
  double det(1.0), num1, num2, total(1.0);

  double* temp = new double[a.rows_];

  for (int i = 0; i < a.rows_; i++) {
    index = i;
    while (index < a.rows_ && a(index, i) == 0) index++;
    if (index == a.rows_) continue;
    if (index != i) {
      for (int j = 0; j < a.rows_; j++) std::swap(a(index, j), a(i, j));

      det = det * pow(-1, index - i);
    }
    for (int j = 0; j < a.rows_; j++) temp[j] = a(i, j);
    for (int j = i + 1; j < a.rows_; j++) {
      num1 = temp[i];
      num2 = a(j, i);
      for (int k = 0; k < a.rows_; k++)
        a(j, k) = (num1 * a(j, k)) - (num2 * temp[k]);
      total = total * num1;
    }
  }

  for (int i = 0; i < a.rows_; i++) det = det * a(i, i);

  delete[] temp;

  return (det / total);
}

S21Matrix S21Matrix::InverseMatrix() const {
  if (Determinant() == 0)
    throw std::invalid_argument("the determinant of the matrix is 0");

  double detKoef = 1.0 / Determinant();
  S21Matrix res = CalcComplements();
  res = res.Transpose();

  return res *= detKoef;
}

S21Matrix& S21Matrix::operator=(const S21Matrix& matrix) {
  if (this == &matrix) return *this;

  Clear();
  rows_ = matrix.rows_;
  cols_ = matrix.cols_;
  Create();
  Copy(matrix);

  return *this;
}

S21Matrix& S21Matrix::operator=(S21Matrix&& matrix) noexcept {
  if (this == &matrix) return *this;

  Clear();
  Move(std::move(matrix));

  return *this;
}

bool S21Matrix::operator==(const S21Matrix& matrix) const {
  return EqMatrix(matrix);
}

S21Matrix& S21Matrix::operator+=(const S21Matrix& matrix) {
  return SumMatrix(matrix);
}

S21Matrix& S21Matrix::operator-=(const S21Matrix& matrix) {
  return SubMatrix(matrix);
}

S21Matrix& S21Matrix::operator*=(const S21Matrix& matrix) {
  return MulMatrix(matrix);
}

S21Matrix& S21Matrix::operator*=(const double num) { return MulNumber(num); }

double& S21Matrix::operator()(const int i, const int j) {
  ThrowOutOfRange(i, j);

  return matrix_[i][j];
}

const double& S21Matrix::operator()(const int i, const int j) const {
  ThrowOutOfRange(i, j);

  return matrix_[i][j];
}

S21Matrix operator+(const S21Matrix& a, const S21Matrix& b) {
  S21Matrix res = a;
  res.SumMatrix(b);

  return res;
}

S21Matrix operator-(const S21Matrix& a, const S21Matrix& b) {
  S21Matrix res = a;
  res.SubMatrix(b);

  return res;
}

S21Matrix operator*(const S21Matrix& a, const S21Matrix& b) {
  S21Matrix res = a;
  res.MulMatrix(b);

  return res;
}

S21Matrix operator*(const S21Matrix& matrix, const double num) {
  S21Matrix res = matrix;
  res.MulNumber(num);

  return res;
}

S21Matrix operator*(const double num, const S21Matrix& matrix) {
  return matrix * num;
}

/*std::ostream& operator<<(std::ostream& out, const S21Matrix& matrix) {
  for (int i = 0; i < matrix.rows_; ++i) {
    for (int j = 0; j < matrix.cols_; ++j) out << matrix(i, j) << " ";
    out << "\n";
  }
  return out;
}*/

void S21Matrix::ThrowDifferentDimensions(const S21Matrix& matrix) const {
  if ((rows_ != matrix.rows_) || (cols_ != matrix.cols_))
    throw std::invalid_argument("different dimensions of matrices");
}

void S21Matrix::ThrowIsNotSquare() const {
  if (rows_ != cols_) throw std::invalid_argument("matrix is not square");
}

void S21Matrix::ThrowOutOfRange(const int i, const int j) const {
  if (!((i >= 0 && i <= rows_ - 1) && (j >= 0 && j <= cols_ - 1)))
    throw std::out_of_range("accessing an element out of range");
}
