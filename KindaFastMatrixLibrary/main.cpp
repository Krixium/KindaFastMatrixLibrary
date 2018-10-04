#include <iostream>
#include <memory>

#include "KindaFastMatrix.h"

#define TEST_N 2

void test()
{
	double *aData = new double[TEST_N * TEST_N];
	double *bData = new double[TEST_N * TEST_N];

	for (int i = 0; i < TEST_N * TEST_N; i++)
	{
		aData[i] = i + 1;
		bData[i] = i + 1;
	}

	std::unique_ptr<kfml::KindaFastMatrix> A = std::make_unique<kfml::KindaFastMatrix>(aData, TEST_N, TEST_N);
	std::unique_ptr<kfml::KindaFastMatrix> B = std::make_unique<kfml::KindaFastMatrix>(bData, TEST_N, TEST_N);
	std::unique_ptr<kfml::KindaFastMatrix> I = std::make_unique<kfml::KindaFastMatrix>(TEST_N);

	std::unique_ptr<kfml::KindaFastMatrix> C;
	C.reset(A->CrossMultiply(B));

	double det = A->GetDeterminant();

	C->Print();
	C->Scale(2);
	C->Print();

	std::cin.get();
}

int main(int argc, char *argv[])
{
	test();

	return 0;
}