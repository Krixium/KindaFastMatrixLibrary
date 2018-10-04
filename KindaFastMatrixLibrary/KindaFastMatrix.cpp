#include "KindaFastMatrix.h"

kfml::KindaFastMatrix::KindaFastMatrix(const size_t m)
	: N(m)
	, M(m)
{
	mData.reset(new double[m * m]);

	for (size_t i = 0; i < m; i++)
	{
		for (size_t j = 0; j < m; j++)
		{
			if (i == j) SetVal(1, i, j);
			else SetVal(0, i, j);
		}
	}
}

kfml::KindaFastMatrix::KindaFastMatrix(double *data, const size_t m, const size_t n)
	: N(n)
	, M(m)
	, mData(data)
	, mInverse(nullptr)
{
}

kfml::KindaFastMatrix *kfml::KindaFastMatrix::CrossMultiply(const KindaFastMatrix& b) const
{
	assert(this->N == b.M);
	size_t K = N;

	KindaFastMatrix *c = new KindaFastMatrix(new double[M * b.N], M, b.N);

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

kfml::KindaFastMatrix *kfml::KindaFastMatrix::CrossMultiply(const std::unique_ptr<KindaFastMatrix>& b) const
{
	return CrossMultiply(*b);
}

void kfml::KindaFastMatrix::Scale(const double scalar)
{
	for (size_t i = 0; i < M * N; i++)
		mData[i] *= scalar;
}

const kfml::KindaFastMatrix *kfml::KindaFastMatrix::GetInverse()
{
	if (mbDetCalculated)
	{
		if (mDet == 0)
		{
			return nullptr;
		}
	}

	if (mInverse == nullptr)
	{
		calculateInverse();
	}

	return mInverse.get();
}

const double kfml::KindaFastMatrix::GetDeterminant()
{
	if (!mbDetCalculated)
	{
		mDet = calculateDet(this, M);
	}

	return mDet;
}

void kfml::KindaFastMatrix::calculateInverse()
{
}

double kfml::KindaFastMatrix::calculateDet(const kfml::KindaFastMatrix *matrix, const size_t size)
{
	double det = 0;
	
	if (matrix->M == 1) return matrix->GetVal(0, 0);

	int sign = 1;
	kfml::KindaFastMatrix tmp(size);
	
	for (size_t i = 0; i < size; i++)
	{
		size_t a = 0;
		size_t b = 0;

		for (size_t row = 0; row < matrix->M; row++)
		{
			for (size_t col = 0; col < matrix->M; col++)
			{
				if (row != 0 && col != i)
				{
					tmp.SetVal(matrix->GetVal(row, col), a, b++);
					if (b == matrix->M - 1)
					{
						a = 0;
						b++;
					}
				}
			}
		}

		det += sign * matrix->GetVal(0, i) * calculateDet(&tmp, size - 1);
		
		sign = -sign;
	}

	return det;
}