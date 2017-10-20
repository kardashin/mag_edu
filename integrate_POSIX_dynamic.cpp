#define HAVE_STRUCT_TIMESPEC T
#include <iostream>
#include <pthread.h>
#include <fstream>
#include <stdlib.h>
#include <ctime>

using namespace std;

pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;

int cycleCounter = 0;	// глобальный счётчик циклов

struct integral_parameters {
	double threadTrapeziumWidth;	// ширина трапеции
	double lowerIntegrationLimit;	// точка начала вычисления интеграла
	int chunkSize;	// число итераций, которые будет проводить каждый поток при каждом запуске цикла while
	int threshold;	// число пробегов цикла while
};

double my_exp(double x) {	// для более долгого вычисления экспоненты :)
	double res = 1., term = 1.;
	int nterms = (x > 0.25 && x < 0.75) ? 1500 : 500;

	for (int n = 1; n <= nterms; n++) {
		term *= x / n;
		res += term;
	}

	return res;
}

double f(double x) {
	return x*my_exp(x);
}

void* integrate(void *recievedData) {	// функция вычисления интеграла; аргумент - указатель на структуру itnerval
	integral_parameters receiveParameters = *(integral_parameters*)recievedData;

	double trapeziumWidth = receiveParameters.threadTrapeziumWidth;
	double lowerLimit = receiveParameters.lowerIntegrationLimit;
	int chunkSize = receiveParameters.chunkSize;
	int threshold = receiveParameters.threshold;
	double result = 0.;	// результат интегрирования

	while (true) {
		pthread_mutex_lock(&mutex);	// для исключения проблемы с вычислением последнего куска интеграла
		++cycleCounter;	// каждый поток увеличивает счётчик на единицу
		if (cycleCounter < threshold + 1) {
			int startStep = (cycleCounter - 1) * chunkSize;	// начальный шаг интегрирования для потока
			pthread_mutex_unlock(&mutex);
			double startPoint = lowerLimit + startStep * trapeziumWidth;	// начальная точка интегрирования для потока
			double functionValueSum = 0.;	// сумма значений функции f()
			int j;
			for (j = 1; j < chunkSize; ++j)	// суммирование значений функций на заданом отрезке за исключением крайних точек
				functionValueSum += f(startPoint + j * trapeziumWidth);

			functionValueSum = (functionValueSum + 0.5 * (f(startPoint) + f(startPoint + j * trapeziumWidth))) * trapeziumWidth;	// применение метода трапеций
			result += functionValueSum;
		}
		else {
			pthread_mutex_unlock(&mutex);
			break;
		}
	}

	return new double(result);
}

int main() {
	const int threadsNumber = 4;	// число потоков
	double lowerLimit = 0., upperLimit = 1.;	// пределы интегрирования
	int defaultStepsNumber = 1000000;	// число точек разбиения интеграла
	int chunkSize = 1000;	// размер блока (число точек), который вычисляет каждый поток

	int newStepsNumber = defaultStepsNumber + chunkSize - defaultStepsNumber % chunkSize;	// переопределение числа разбиений отрезка интегрирования на случай, если оно не делится на число потоков нацело
	double trapeziumWidth = (upperLimit - lowerLimit) / newStepsNumber;	// ширина трапеции, которая будет передана каждому потоку для вычисления интеграла в точке
	int threshold = newStepsNumber / chunkSize;	// число итераций на один поток

	pthread_t threadID[threadsNumber];	// структура, содержащая параметры интегрирования для каждого потока
	integral_parameters sendParameters[threadsNumber];

	for (int i = 0; i < threadsNumber; i++) {	// передача аргументов функции integrate() через структуру sendParameters
		sendParameters[i].threadTrapeziumWidth = trapeziumWidth;
		sendParameters[i].lowerIntegrationLimit = lowerLimit;
		sendParameters[i].chunkSize = chunkSize;
		sendParameters[i].threshold = threshold;

		pthread_create(&threadID[i], 0, integrate, &sendParameters[i]);
	}

	cout.precision(15);
	double finalResuilt = 0.;	// результат интегрирования
	for (int i = 0; i < threadsNumber; i++) {
		double* threadResult;
		pthread_join(threadID[i], (void**)&threadResult);
		cout << "Thread number " << i + 1 << " integration result: " << *threadResult << endl; //" complted." << endl;
		finalResuilt += *threadResult;	// суммирование результатов интегрирования всех потоков
		delete threadResult;
	}

	cout << endl << "Final result: " << finalResuilt << endl << "Execution time: " << clock() << " ms." << endl;
	return 0;
}
