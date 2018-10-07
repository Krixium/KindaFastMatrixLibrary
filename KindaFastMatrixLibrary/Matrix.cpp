#include "Matrix.h"

kfml::Matrix::Matrix(const size_t m)
	: M(m)
	, N(m)
	, mData(std::make_unique<std::vector<double>>(m * m))
	, mInverse(nullptr)
	, mDet(0)
	, mbDetCalculated(false)
{
	for (size_t i = 0; i < m; i++)
		SetVal(1, i, i);
}

kfml::Matrix::Matrix(const size_t m, const size_t n)
	: M(m)
	, N(n)
	, mData(std::make_unique<std::vector<double>>(m * n))
	, mInverse(nullptr)
	, mDet(0)
	, mbDetCalculated(false)
{
	for (size_t i = 0; i < m; i++)
		for (size_t j = 0; j < n; j++)
			SetVal(1, i, j);
}

kfml::Matrix::Matrix(double *data, const size_t m, const size_t n)
	: N(n)
	, M(m)
	, mData(std::make_unique<std::vector<double>>())
	, mInverse(nullptr)
	, mDet(0)
	, mbDetCalculated(false)
{
	for (size_t i = 0; i < m * n; i++)
		mData->push_back(data[i]);
}

kfml::Matrix *kfml::Matrix::CrossMultiply(const Matrix& b)
{
	assert(N == b.M);

	size_t K = N;

	Matrix *c = new Matrix(new double[M * b.N], M, b.N);

	for (size_t i = 0; i < M; i++)
	{
		for (size_t j = 0; j < b.N; j++)
		{
			double sum = 0;
			for (size_t k = 0; k < K; k++)
			{
				sum += GetVal(i, k) * b.GetVal(k, j);
			}
			c->SetVal(sum, i, j);
		}
	}

	return c;
}

void kfml::Matrix::Scale(const double scalar)
{
	for (size_t i = 0; i < M * N; i++)
		mData->at(i) *= scalar;
}

kfml::Matrix *kfml::Matrix::GetInverse()
{
	if (!mbDetCalculated)
	{
		GetDeterminant();
	}

	if (mDet == 0)
	{
		return nullptr;
	}

	if (mInverse == nullptr)
	{
		calculateInverse();
	}

	return mInverse.get();
}

const double kfml::Matrix::GetDeterminant()
{
	if (!mbDetCalculated)
	{
		mDet = calculateDet(*this, M);
		mbDetCalculated = true;
	}

	return mDet;
}

void kfml::Matrix::calculateInverse()
{
	assert(M == N);

	mInverse = std::make_unique<kfml::Matrix>(M, N);

	kfml::Matrix tmp(M - 1, N - 1);

	double inverseDet = 1 / GetDeterminant();
	int sign;

	for (size_t i = 0; i < M; i++)
	{
		for (size_t j = 0; j < N; j++)
		{
			(i + j + 2) % 2 == 0 ? sign = 1 : sign = -1;
			extractMinor(*this, tmp, i, j, M);
			mInverse->SetVal(inverseDet * sign * tmp.GetDeterminant(), i, j);
		}
	}

	mInverse->Transpose();
}

double kfml::Matrix::calculateDet(const kfml::Matrix &matrix, const size_t size)
{
	double det = 0;
	
	if (matrix.M == 1) return matrix.GetVal(0, 0);

	int sign = 1;
	kfml::Matrix tmp(size - 1);
	
	for (size_t i = 0; i < size; i++)
	{
		extractMinor(matrix, tmp, 0, i, size);
		det += sign * matrix.GetVal(0, i) * calculateDet(tmp, size - 1);
		sign = -sign;
	}

	return det;
}

void kfml::Matrix::extractMinor(const kfml::Matrix& matrix, kfml::Matrix& tmp, const size_t x, const size_t y, const size_t n)
{
	assert(matrix.M == matrix.N);
	assert(x <= matrix.M - 1);
	assert(y <= matrix.N - 1);
	assert(n <= matrix.M && n <= matrix.N);

	int i = 0;
	int j = 0;

	for (size_t row = 0; row < n; row++)
	{
		for (size_t col = 0; col < n; col++)
		{
			if (row != x && col != y)
			{
				tmp.SetVal(matrix.GetVal(row, col), i, j++);
				if (j == matrix.M - 1)
				{
					j = 0;
					i++;
				}
			}
		}
	}
}
