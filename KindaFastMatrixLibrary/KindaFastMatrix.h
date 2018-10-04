#pragma once

#include <assert.h>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace kfml
{
	class KindaFastMatrix
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
		std::unique_ptr<double[]> mData; 

		// Determinant
		double mDet; 
		// Inverse Matrix
		std::unique_ptr<KindaFastMatrix> mInverse; 

	public:
		KindaFastMatrix(const size_t m);
		KindaFastMatrix(double *data, const size_t m, const size_t n);

		KindaFastMatrix *CrossMultiply(const KindaFastMatrix& b) const;
		KindaFastMatrix *CrossMultiply(const std::unique_ptr<KindaFastMatrix>& b) const;
		void Scale(const double scalar);
		const KindaFastMatrix *GetInverse();
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

		inline double GetVal(const size_t m, const size_t n) const
		{
			return mData[m * N + n];
		}

		inline void SetVal(const double val, const size_t m, const size_t n)
		{
			mData[m * N + n] = val;
			mbDetCalculated = false;
		}

		inline void Print()
		{
			for (size_t i = 0; i < M; i++)
			{
				for (size_t j = 0; j < N; j++)
				{
					std::cout << std::to_string((int)GetVal(i, j)) << "\t";
				}
				std::cout << std::endl << std::endl;
			}
			std::cout << std::endl;
		}

	private:
		void calculateInverse();
		static double calculateDet(const KindaFastMatrix *matrix, const size_t size);
		static void extractCofactor(const KindaFastMatrix& matrix, KindaFastMatrix& tmp, const size_t x, const size_t y, const size_t n);

		// OpenCL helpers here
	};

}
