CXX = g++
CXXFLAGS = -Wall -Werror -Wextra -pedantic -std=c++17
LDFLAGS = -lgtest --coverage

all: s21_matrix_oop.a

build: s21_matrix_oop.cc
	$(CXX) $(CXXFLAGS) -c s21_matrix_oop.cc -o s21_matrix_oop.o

s21_matrix_oop.a: build
	ar rcs s21_matrix_oop.a s21_matrix_oop.o

clean:
	rm -rf *.o *.out s21_matrix_oop.a test
	rm -rf *.gcda *.gcno report gcov_report.* gcov_report *.info

rebuild: clean all

leaks: test
	leaks -atExit -- ./test

clang:
	cp ../materials/linters/.clang-format ../src/.clang-format
	clang-format -style=Google -n *.cc *.h
	rm .clang-format

test: s21_matrix_oop.a ./test.cc
	$(CXX) $(CXXFLAGS) $(LDFLAGS) s21_matrix_oop.a s21_matrix_oop.cc ./test.cc -o test
	./test

gcov_report: clean test
	lcov -no-external -t "stest" -o s21_test.info -c -d .
	genhtml -o report s21_test.info
	open ./report/index.html
