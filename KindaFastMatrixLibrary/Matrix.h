#pragma once

#include <assert.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace kfml
{
	class Matrix
	{
	public:
		// Row
		const size_t M;	
		// Column
		const size_t N;	

	private:
		// Was the determinant already calculated?
		bool mbDetCalculated; 

		// Matrix data
		std::unique_ptr<std::vector<double>> mData;

		// Determinant
		double mDet; 
		// Inverse Matrix
		std::unique_ptr<Matrix> mInverse;

	public:
		Matrix(const size_t m);
		Matrix(const size_t m, const size_t n);
		Matrix(double *data, const size_t m, const size_t n);

		Matrix *CrossMultiply(const Matrix& b);
		inline Matrix *CrossMultiply(const std::unique_ptr<Matrix>& b) { return CrossMultiply(*b); }
		inline Matrix *CrossMultiply(const Matrix *b) { return CrossMultiply(*b); }
		void Scale(const double scalar);
		Matrix *GetInverse();
		const double GetDeterminant();

		inline void Transpose()
		{
			for (size_t i = 0; i < M; i++)
			{
				for (size_t j = i + 1; j < N; j++)
				{
					double tmp = GetVal(i, j);
					SetVal(GetVal(j, i), i, j);
					SetVal(tmp, j, i);
				}
			}
		}

		inline double GetVal(const size_t m, const size_t n) const { return mData->at(m * N + n); }

		inline void SetVal(const double val, const size_t m, const size_t n)
		{
			mData->at(m * N + n) = val;
			mbDetCalculated = false;
			mInverse.reset();
		}

		inline void Print()
		{
			for (size_t i = 0; i < M; i++)
			{
				for (size_t j = 0; j < N; j++)
				{
					std::cout << std::to_string(GetVal(i, j)) << "\t";
				}
				std::cout << std::endl << std::endl;
			}
			std::cout << std::endl;
		}

		inline void PrintLine()
		{
			std::cout << "[";
			for (size_t i = 0; i < M * N - 1; i++)
			{
				std::cout << std::round(mData->at(i)) << ",";
			}
			std::cout << std::round(mData->at(M * N - 1)) << "]";
			std::cout << std::endl;
		}

		inline void ZeroOut() { for (size_t i = 0; i < M * N; i++) mData->at(i) = 0; }

	private:
		void calculateInverse();
		static double calculateDet(const Matrix &matrix, const size_t size);
		static void extractMinor(const Matrix& matrix, Matrix& tmp, const size_t x, const size_t y, const size_t n);
	};
}
