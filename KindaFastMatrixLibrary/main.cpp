#include <algorithm>
#include <chrono>
#include <functional>
#include <fstream>
#include <iostream>
#include <iterator>
#include <memory>
#include <string>
#include <sstream>
#include <utility>
#include <vector>

#include "Matrix.h"

const size_t LINE_LIMIT = 8096;

// Time things in seconds 
double TimeExecution(std::function<void()> callback)
{
	std::chrono::time_point<std::chrono::steady_clock> mStartTime = std::chrono::high_resolution_clock::now();
	callback();
	std::chrono::time_point<std::chrono::steady_clock> mEndTime = std::chrono::high_resolution_clock::now();

	std::chrono::duration<double> delta = mEndTime - mStartTime;
	return delta.count();
}

// Extracts the size of the matrix from the file
size_t parseMatrixSize(char *lineBuffer, const size_t limit)
{
	for (size_t i = 2; i < limit; i++)
	{
		if (lineBuffer[i] == 'x')
		{
			lineBuffer[i] = 0;
			break;
		}
	}

	return atoi(lineBuffer + 1);
}

// Splits a string
void split(const std::string& str, std::vector<double>& cont, char delim = ' ')
{
	std::stringstream ss(str);
	std::string token;
	while (std::getline(ss, token, delim)) {
		cont.push_back(stod(token));
	}
}

// Parses the file into pairs of matrix A and matrix B and shoves it all into a vector
std::vector<std::pair<kfml::Matrix, kfml::Matrix>> *parseFile(const std::string& filename)
{
	std::ifstream dataFile(filename);

	if (!dataFile)
	{
		return nullptr;
	}

	std::vector<std::pair<kfml::Matrix, kfml::Matrix>> *data = new std::vector<std::pair<kfml::Matrix, kfml::Matrix>>();

	size_t matrixSize = 0;
	char lineBuffer[LINE_LIMIT];

	while (dataFile)
	{
		dataFile.getline(lineBuffer, LINE_LIMIT);

		if (dataFile.eof())
		{
			continue;
		}
		if (lineBuffer[0] == '#')
		{
			matrixSize = parseMatrixSize(lineBuffer, LINE_LIMIT);
		}
		else
		{
			// Grab all the values of the matrix and put into a single vector
			if (matrixSize == 0) continue;

			std::vector<double> cont;
			split(lineBuffer, cont, ',');

			for (size_t i = 0; i < matrixSize - 1; i++)
			{
				dataFile.getline(lineBuffer, LINE_LIMIT);
				split(lineBuffer, cont, ',');
			}

			// AX = B
			double *A = new double[matrixSize * matrixSize];
			double *B = new double[matrixSize];
			size_t indexA = 0;
			size_t indexB = 0;

			// Split into A and B
			for (size_t i = 0; i < cont.size(); i++)
			{
				if (i % (matrixSize + 1) == matrixSize) B[indexB++] = cont[i];
				else A[indexA++] = cont[i];
			}

			data->push_back(std::pair<kfml::Matrix, kfml::Matrix>(kfml::Matrix(A, matrixSize, matrixSize), kfml::Matrix(B, matrixSize, 1)));
			matrixSize = 0;
		}
	}

	dataFile.close();
	return data;
}

int main(int argc, char *argv[])
{
	if (argc != 2)
	{
		std::cout << "Incorrect number of commandline arguements" << std::endl;
		return 0;
	}

	double time;
	std::unique_ptr<std::vector<std::pair<kfml::Matrix, kfml::Matrix>>> inputData;

	time = TimeExecution([&]() {
		inputData.reset(parseFile(argv[1]));
	});
	std::cout << "Reading the file took: " << time << "s" << std::endl;

	if (inputData == nullptr)
	{
		std::cout << "Could not read data in from file";
		return 0;
	}

	for (size_t i = 0; i < inputData->size(); i++)
	{
		kfml::Matrix& A = inputData->at(i).first;
		kfml::Matrix& B = inputData->at(i).second;

		time = TimeExecution([&]() {
			A.GetInverse();
		});
		std::cout << "Calculating inverse took: " << time << "s" << std::endl;

		std::unique_ptr<kfml::Matrix> C;

		time = TimeExecution([&]() {
			C.reset(A.GetInverse()->CrossMultiply(B));
		});
		std::cout << "Solving the system of equations took: " << time << "s" << std::endl;

		C->PrintLine();
	}

	//std::cout << "Press any key to exit" << std::endl;
	//std::cin.get();
	return 0;
}