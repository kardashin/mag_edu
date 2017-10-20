#include <cmath>
#include <iostream>
#include <mpi.h>
#include <ctime>
#include <fstream>

using namespace std;	

double my_exp(double x) {	// для более долгого вычисления экспоненты :)
	double res = 1., term = 1.;
	int nterms = 1000;

	for (int n = 1; n <= nterms; n++) {
		term *= x / n;
		res += term;
	}

	return res;
}

double f(double x) {
	return x*my_exp(x);
}

double integrate(int startStep, int chunkSize, double lowerLimit, double trapeziumWidth) {	// функция интегрирования
	double result = 0.;	// результат интегрирования
	double startPoint = lowerLimit + startStep * trapeziumWidth;	// точка начала интегрирования
	
	int i;
	for (i = 1; i < chunkSize; ++i)
		result += f(startPoint + i * trapeziumWidth);	// суммирование значений f() в точках отрезка, исключая начальную и конечную
	
	result = (result + 0.5 * (f(startPoint) + f(startPoint + i * trapeziumWidth))) * trapeziumWidth;	// применение метода трапеций
	
	return result;
}

int main(int argc, char** argv) {
	double lowerLimit = 0., upperLimit = 1.;	// верхний и нижний пределы интегрирования
	int stepsNumber = 1000000;	// число разбиений отрезка интегрирования
	int chunkSize = 1000;	// число подзадач на каждый поток
	double threadIntegrationResult = 0.;
	int cycleCounter = 0;	// счётчик циклов
	
	double trapeziumWidth = (upperLimit - lowerLimit) / stepsNumber;	// ширина трапеции, которая будет передана каждому потоку для вычисления интеграла в точке
	int threshold = stepsNumber / chunkSize;	// число итераций на один поток
	
	if(MPI_Init(&argc, &argv))
		return -1;
    int size, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Status status;
	
	bool threadResponse;	// переменная для фиксации факта ответа вычисляющих потоков на запросы главного потока
	int completionFlag = 0;	// флаг завершения вычисления
	if (rank == 0) {	// главный поток, раздающий задачи другим процессам
		for (cycleCounter; cycleCounter < threshold; ++cycleCounter) {	// распределение задач потокам
			MPI_Recv(&threadResponse, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);	// ожидание ответа от потока
			MPI_Send(&cycleCounter, 1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD);	// посылка счетчика циклов
		}

		completionFlag = -1;	// изменение значения флага; передпча флага начинает завершение работы вычисляющих процессов
		
		for (int i = 0; i < size - 1; ++i) {	// завершение работы вычисляющих потоков
			MPI_Recv(&threadResponse, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			MPI_Send(&completionFlag, 1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD);	// передача флага завершения вычислений
		}
	}
	else	// потоки, выполняющие интегрирование
		while (cycleCounter != -1) {	// работает пока не будет отправлен completionFlag = -1
			MPI_Send(&threadResponse, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);	// посылка ответа главному процессу
			MPI_Recv(&cycleCounter, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);	// получение от главного процесса счётчика циклов
			if (cycleCounter != -1)
				threadIntegrationResult += integrate(cycleCounter * chunkSize, chunkSize, lowerLimit, trapeziumWidth);
		}

	cout.precision(15);
	cout << "Process: " << rank << ". Result: " << threadIntegrationResult << "." << endl;
		
	double* gatheredData = new double [size];   // выделение памяти под массив для сбора данных с потоков
	MPI_Gather(&threadIntegrationResult, 1, MPI_DOUBLE, gatheredData, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);	// сбор данных со всех процессов

	double finalIntegrationResult = 0.;
	for (int i = 0; i < size; ++i)
		finalIntegrationResult += gatheredData[i];

	MPI_Finalize();

	if (rank == 0)
		cout << endl << "Final result: " << finalIntegrationResult << "." << endl;
	
	return 0;
}
