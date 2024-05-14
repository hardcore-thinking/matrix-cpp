#include <iostream>
#include <cstdlib>
#include <vector>
#include <random>

double GetRandomDouble(double a = 0.0, double b = 1.0) {
	static std::random_device randomDevice;
	std::uniform_real_distribution<double> distribution(a, b);
	return distribution(randomDevice);
}

template <size_t N, size_t M>
class Matrix;

template <size_t N, size_t M>
std::ostream& operator<<(std::ostream& out, Matrix<N, M> const& mat);

template <size_t N, size_t M>
Matrix<N, M> operator+(Matrix<N, M> const& lhs, Matrix<N, M> const& rhs);

template <size_t N, size_t M>
Matrix<N, M> operator-(Matrix<N, M> const& mat);

template <size_t N, size_t M>
Matrix<N, M> operator-(Matrix<N, M> const& lhs, Matrix<N, M> const& rhs);

template <size_t N, size_t M>
Matrix<N, M> operator*(double k, Matrix<N, M> const& mat);

template <size_t N, size_t M>
Matrix<N, M> operator*(Matrix<N, M> const& mat, double k);

template <size_t N, size_t M = N>
class Matrix {
	public:
		Matrix();
		Matrix(Matrix const&);

		size_t Dimension() const;
		size_t Lines() const;
		size_t Columns() const;

		static Matrix Random();
		static Matrix Identity();
		static Matrix<M, N> GetTranspose(Matrix<N, M>& mat);

		void operator=(Matrix const& mat);
		
		template <size_t M2, size_t N2> friend class Matrix;

		friend std::ostream& operator<<<>(std::ostream& out, Matrix const& mat);
		friend Matrix operator+<>(Matrix const& lhs, Matrix const& rhs);
		friend Matrix operator-<>(Matrix const& mat);
		friend Matrix operator-<>(Matrix const& lhs, Matrix const& rhs);
		friend Matrix operator*<>(double k, Matrix const& mat);
		friend Matrix operator*<>(Matrix const& mat, double k);
		
	private:
		std::vector<double> _elements = {}; 
		size_t _n = 0;
		size_t _m = 0;
};

template <size_t N, size_t M>
Matrix<N, M>::Matrix() : _n(N), _m(M) {
	_elements.resize(_n * _m);
}

template <size_t N, size_t M>
Matrix<N, M>::Matrix(Matrix const& mat) : _n(mat._n), _m(mat._m) {
	_elements.resize(mat._n, mat._m);
	for (size_t i = 0; i < mat._n; i++) {
		for (size_t j = 0; j < mat._m; j++) {
			_elements[i * M + j] = mat._elements[i * M + j];
		}
	}
}

template <size_t N, size_t M>
size_t Matrix<N, M>::Dimension() const {
	return _n * _m;
}

template <size_t N, size_t M>
size_t Matrix<N, M>::Lines() const {
	return _n;
}

template <size_t N, size_t M>
size_t Matrix<N, M>::Columns() const {
	return _m;
}

template <size_t N, size_t M>
Matrix<N, M> Matrix<N, M>::Random() {
	Matrix<N, M> result;

	for (size_t i = 0; i < N; i++) {
		for (size_t j = 0; j < M; j++) {
			result._elements[i * M + j] = std::round(GetRandomDouble(0, 9));
		}
	}

	return result;
}

template <size_t N, size_t M>
Matrix<N, M> Matrix<N, M>::Identity() {
	static_assert(N == M, "An identity matrix only exists for square matrix.");
	
	Matrix<N> identity;

	for (size_t i = 0; i < N; i++) {
		identity._elements[N * i + i] = 1;
	}

	return identity;
}

template <size_t N, size_t M>
Matrix<M, N> Matrix<N, M>::GetTranspose(Matrix<N, M>& mat) {
	Matrix<M, N> transposed;

	for (size_t i = 0; i < N; i++) {
		for (size_t j = 0; j < M; j++) {
			transposed._elements[j * N + i] = mat._elements[i * M + j];
		}
	}

	return transposed;
}

template <size_t N, size_t M>
void Matrix<N, M>::operator=(Matrix<N, M> const& mat) {
	_n = mat._n;
	_m = mat._m;
	
	for (size_t i = 0; i < mat._n; i++) {
		for (size_t j = 0; j < mat._m; j++) {
			_elements[i * M + j] = mat._elements[i * M + j];
		}
	}

}

template <size_t N, size_t M>
std::ostream& operator<<(std::ostream& out, Matrix<N, M> const& mat) {
	for (size_t i = 0; i < mat._n; i++) {
		if (i != 0) {
			std::cout << std::endl;
		}

		for (size_t j = 0; j < mat._m; j++) {
			out << mat._elements[i * mat._m + j] << " ";
		}
	}

	return out;
}

template <size_t N, size_t M>
Matrix<N, M> operator+(Matrix<N, M> const& lhs, Matrix<N, M> const& rhs) {
	Matrix<N, M> result;

	for (size_t i = 0; i < N; i++) {
		for (size_t j = 0; j < M; j++) {
			result._elements[i * M + j] = lhs._elements[i * M + j] + rhs._elements[i * M + j];
		}
	}

	return result;
}

template <size_t N, size_t M>
Matrix<N, M> operator-(Matrix<N, M> const& mat) {
	Matrix<N, M> result;

	for (size_t i = 0; i < N; i++) {
		for (size_t j = 0; j < M; j++) {
			result._elements[i * M + j] = (-1) * mat._elements[i * M + j];
		}
	}

	return result;
}

template <size_t N, size_t M>
Matrix<N, M> operator-(Matrix<N, M> const& lhs, Matrix<N, M> const& rhs) {
	return lhs + (-rhs);
}

template <size_t N, size_t M>
Matrix<N, M> operator*(double k, Matrix<N, M> const& mat) {
	Matrix<N, M> result;

	for (size_t i = 0; i < N; i++) {
		for (size_t j = 0; j < M; j++) {
			result._elements[i * M + j] = k * mat._elements[i * M + j];
		}
	}

	return result;

}

template <size_t N, size_t M>
Matrix<N, M> operator*(Matrix<N, M> const& mat, double k) {
	Matrix<N, M> result = k * mat;

	return result;
}

int main() {
	Matrix<2> matrix22;

	Matrix<3> matrix33 = Matrix<3>::Identity();

	Matrix<3, 4> matrix34 = Matrix<3, 4>::Random();

	std::cout << "Dimension: " << matrix22.Dimension() << std::endl;
	std::cout << "Lines:     " << matrix22.Lines()     << " lines"   << std::endl;
	std::cout << "Columns:   " << matrix22.Columns()   << " columns" << std::endl << std::endl;

	std::cout << matrix33 << std::endl << std::endl;

	std::cout << Matrix<3>::Identity() << std::endl << std::endl;

	std::cout << matrix34 << std::endl << std::endl;

	std::cout << Matrix<3, 4>::GetTranspose(matrix34) << std::endl << std::endl;

	std::cout << Matrix<2>::Identity() + Matrix<2>::Identity() << std::endl << std::endl;

	std::cout << -Matrix<4>::Identity() << std::endl << std::endl;

	std::cout << Matrix<2>::Identity() - Matrix<2>::Identity() << std::endl << std::endl;

	std::cout << 3 * Matrix<3>::Identity() << std::endl << std::endl;

	std::cout << Matrix<3>::Identity() * 3 << std::endl << std::endl;

	return EXIT_SUCCESS;
}
