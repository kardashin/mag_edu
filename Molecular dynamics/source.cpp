#define _USE_MATH_DEFINES	// использование математических констант
#include <iostream>
#include <fstream>
#include <cstdio>
#include <vector>
#include <cmath>
#include <random>

using namespace std;

// хороший генератор случайных чисел от 0 до 1
random_device                  rand_dev;
mt19937                        generator(rand_dev());
uniform_int_distribution<int>  distr(0, 1);

class VectorComponents {	// класс векторов в трёхмерном пространстве; перегрузка операторов для действий над векторами
public:
	double x;
	double y;
	double z;

	VectorComponents operator+(const VectorComponents& summand) // сложение векторов
	{
		VectorComponents temp;
		temp.x = x + summand.x;
		temp.y = y + summand.y;
		temp.z = z + summand.z;
		return temp;
	}

	VectorComponents operator-(const VectorComponents& subtrahend)	// вычитание векторов
	{
		VectorComponents temp;
		temp.x = x - subtrahend.x;
		temp.y = y - subtrahend.y;
		temp.z = z - subtrahend.z;
		return temp;
	}

	VectorComponents operator*(double multiplier)	// умножение вектора на число
	{
		VectorComponents temp;
		temp.x = x * multiplier;
		temp.y = y * multiplier;
		temp.z = z * multiplier;
		return temp;
	}

	VectorComponents operator/(double devider)	// деление вектора на число
	{
		VectorComponents temp;
		temp.x = x / devider;
		temp.y = y / devider;
		temp.z = z / devider;
		return temp;
	}

	VectorComponents operator=(double assignment)	// присвоение компонент
	{
		x = assignment;
		y = assignment;
		z = assignment;
		return *this;
	}

	VectorComponents& operator+=(const VectorComponents& summand) {
		x = x + summand.x;
		y = y + summand.y;
		z = z + summand.z;
		return *this;
	}

	VectorComponents& operator-=(const VectorComponents& subtrahend) {
		x = x - subtrahend.x;
		y = y - subtrahend.y;
		z = z - subtrahend.z;
		return *this;
	}

	VectorComponents operator%(double denominator)
	{
		VectorComponents temp;
		temp.x = remainder(x, denominator);
		temp.y = remainder(y, denominator);
		temp.z = remainder(z, denominator);
		return temp;
	}

	double SquaredVectorNorm() {	// вычисление квадрата нормы вектора
		return x*x + y*y + z*z;
	}
};

const double boltzmannСonstant = 1;	// постоянная Больцмана
vector <VectorComponents> position;	// вектор радиус-векторов частиц
vector <VectorComponents> velocity;	// вектор скоростей частиц
vector <VectorComponents> force;	// вектор сил, действующих на частицы
vector <int> intervals;	// вектор с количеством интервалов скоростей в распределении Максвелла
double temperature, density, sigma, epsilon, particleMass, timeStep, calculationTime, cutRadius;	// температура, плотность, параметры сигма и эпсилон, масса частиц, шаг по времени, радиус отсечки потенциала
double sigmaPow6, sigmaPow12;	// параметр сигма в 6-й и 12-й степенях
double potentialEnergy = 0, kineticEnergy = 0;	// кинетическая и потенциальная энергии системы
double squaredCutRadius;	// квардрат радиуса отсечки потенциала
double currentTime;	// время, прошедшее с начала вычислений
double lambda;	// параметр лямбда для термостата Бернедсена
double relaxationTime = 10;	// приблизительное время вывода системы в равновесие
double speedInterval;	// интервал скоростей для распределения Максвелла
double maxSpeed;	// скорость, до которой будет вычисляться распределение Максвелла
int intervalsNumber;	// число интервалов при вычислении распределения Максвелла
int numberOfParticles;	// число частиц

void Initialization() {	// функция инициализации
	int numberOfParticlesInTheLine;	// число частиц в одной линии ячейки
	double distanceBetweenParticles;	// расстояние между частицами в ячейке
	int currentParticle = 0;

    // считывание параметров
	ifstream finParameters("inputParameters.txt");
	finParameters >> temperature >> density >> sigma >> epsilon >> particleMass >> timeStep >> calculationTime >> cutRadius >> intervalsNumber;
	finParameters.close();

	sigmaPow6 = pow(sigma, 6);
	sigmaPow12 = pow(sigma, 12);
	squaredCutRadius = pow(cutRadius, 2);
	numberOfParticles = int(pow(2. * cutRadius, 3) * density);
//	numberOfParticles = 64;	// для теста; извлекается кубический корень
	numberOfParticlesInTheLine = int(ceil(cbrt(double(numberOfParticles))));
	distanceBetweenParticles = 2. * cutRadius / numberOfParticlesInTheLine;

	maxSpeed = 3. * sqrt(3. * temperature * boltzmannСonstant / particleMass);
	speedInterval = maxSpeed / intervalsNumber;	// вычисление интервала скоростей

    // изменение размеров векторов
	position.resize(numberOfParticles);
	velocity.resize(numberOfParticles);
	force.resize(numberOfParticles);
	intervals.resize(intervalsNumber);
	
/* плохая (?) реализация расстановки частитц
	int iterationCount = 0;	// счётчик итераций
	double halfCrystalLength = crystalLength / 2.;	// половина длины ячейки
	double startXcoordinate = (-1) * halfCrystalLength + halfCrystalLength / crystalDistance;	// координата x первой частицы в линии ячейки
	double startYcoordinate = startXcoordinate;	// координата y первой частицы в линии ячейки
	double startZcoordinate = startXcoordinate;	// координата z первой частицы в линии ячейки
	double endXcoordinate = halfCrystalLength - halfCrystalLength / crystalDistance;	// координата x последней частицы в линии ячейки
	double endYcoordinate = endXcoordinate;	// координата y последней частицы в линии ячейки
	double endZcoordinate = endXcoordinate;	// координата z последней частицы в линии ячейки
	double placingStep = crystalLength / crystalDistance;	// шаг расстановки частиц в ячейки
	for (double a = startXcoordinate; a <= endXcoordinate; a += placingStep)
		for (double b = startYcoordinate; b <= endYcoordinate; b += placingStep)
			for (double c = startZcoordinate; c <= endZcoordinate; c += placingStep) {
				if (iterationCount == numberOfParticles)
					break;
				position[iterationCount].x = a;
				position[iterationCount].y = b;
				position[iterationCount].z = c;
				++iterationCount;
			}
*/
	// расстановка частиц в узлах решётки
	for (int i = 0; i < numberOfParticlesInTheLine; ++i) {
		for (int j = 0; j < numberOfParticlesInTheLine; ++j) {
			for (int k = 0; k < numberOfParticlesInTheLine; ++k) {
				if (currentParticle >= numberOfParticles) {
					break;
				}
				position[numberOfParticlesInTheLine*numberOfParticlesInTheLine*i + numberOfParticlesInTheLine*j + k].x = distanceBetweenParticles*i;
				position[numberOfParticlesInTheLine*numberOfParticlesInTheLine*i + numberOfParticlesInTheLine*j + k].y = distanceBetweenParticles*j;
				position[numberOfParticlesInTheLine*numberOfParticlesInTheLine*i + numberOfParticlesInTheLine*j + k].z = distanceBetweenParticles*k;
				++currentParticle;
			}
		}
	}

    // присвоение частицам случайных компонент начаных скоростей в интервале от 0 до 1
	for (int i = 0; i < numberOfParticles; ++i) {
		velocity[i].x = distr(generator);
		velocity[i].y = distr(generator);
		velocity[i].z = distr(generator);
	}
}

void CalculateEnergyAndForces() {	// функция вычисления сил и энергий
	VectorComponents distanceVector;	// вектор разности радиус-векторов i-й и j-й частиц
	double squaredDistance;	// квадрат расстояния между частицами
	double systemTemperature;	// мгновенная температура системы

    // обнуление энергий
	kineticEnergy = 0;
	potentialEnergy = 0;

	for (int i = 0; i < numberOfParticles; ++i) {
		force[i] = 0;	// обнуление сил
		kineticEnergy += velocity[i].SquaredVectorNorm() * particleMass / 2;	// вычисление полной кинетической энергии

		for (int j = 0; j < i; ++j) {
			distanceVector = position[i] - position[j];	// вычисление вектора расстояния между i-й и j-й частицами
			distanceVector = distanceVector % (2 * cutRadius);	// граничные условия
			squaredDistance = distanceVector.SquaredVectorNorm();	// вычисление квадрата расстояния между частицами

		    // вспомогательные вычисления
			double distancePow6 = pow(squaredDistance, 3);
			double distancePow12 = pow(distancePow6, 2);
			double distancePow8 = distancePow6 * squaredDistance;
			double distancePow14 = distancePow12 * squaredDistance;

			potentialEnergy += 4. * epsilon * (sigmaPow12 / distancePow12 - sigmaPow6 / distancePow6);	// вычисление полной потенциальной энергии

			if (squaredDistance < squaredCutRadius) {	// учёт радиуса отсечки потенциала
				force[i] += distanceVector * (24. * epsilon * (2. * sigmaPow12 / distancePow14 - sigmaPow6 / distancePow8));	// сила, действующая на i-ю частицу со стороны j-й
				force[j] -= distanceVector * (24. * epsilon * (2. * sigmaPow12 / distancePow14 - sigmaPow6 / distancePow8));	// сила, действующая на j-ю частицу со стороны i-й
			}
		}
	}

	systemTemperature = 2. * kineticEnergy / (3. * boltzmannСonstant * numberOfParticles);	// мгновенная температура системы
	lambda = sqrt(1 + (temperature / systemTemperature - 1.) * timeStep / relaxationTime);
}

void SolveMotionEquation() {	// функция решения уравнения движения
	for (int i = 0; i < numberOfParticles; ++i) {
		velocity[i] = (velocity[i] + force[i] * timeStep / particleMass) * lambda;	// вычисление новых скоростей частиц
		position[i] += velocity[i] * timeStep;	// вычисление новых радиус-векторов частиц
	}
}

ofstream foutEnergy("outputEnergy.txt");
//ofstream foutPosition("outputPosition.txt");	\\ вывод в файл координат всех частиц на каждом шаге
ofstream foutMaxwell("outputMaxwell.txt");
ofstream foutTheorMaxwell("outputTheorMaxwell.txt");

void RecordState() {
	foutEnergy << currentTime << "\t" << kineticEnergy/numberOfParticles << "\t" << potentialEnergy / numberOfParticles << "\t" << (potentialEnergy + kineticEnergy) / numberOfParticles << endl;	// запись времени и энергий в файл
/*																																		
	for (int i = 0; i < numberOfParticles; ++i)	// запись координат частиц в файл
		foutPosition << position[i].x << "\t" << position[i].y << "\t" << position[i].z << "\t" << endl;
	foutPosition << endl;*/
}

void CalculateMaxwellDistribution() {	// функция вычисления распределения Максвелла
	double particleSpeed = 0.;
	double speedInterval = maxSpeed / intervalsNumber;

	for (int i = 0; i < numberOfParticles; ++i) {
		particleSpeed = sqrt(velocity[i].SquaredVectorNorm());	// скорость i-й частицы
		++intervals[int(ceil(particleSpeed / speedInterval))];	// интервал скоростей, в который попала i-я частица
	}

/*
\\ плохая версия
	vector <double> speed(numberOfParticles);	//	вектор со скоростями частиц
	for (int i = 0; i < numberOfParticles; ++i)
		speed[i] = sqrt(velocity[i].SquaredVectorNorm());
	for (int i = 0; i < numberOfParticles; ++i)
		for (int j = 0; j < intervalsNumber; ++j)
			if (speed[i] >= speedInterval * j && speed[i] <= speedInterval * (j + 1))
				++intervals[j];
*/
}

double IntegratedMaxwellDistribution(double t, double v) {	// интеграл от функции плотности распределения Максвелла
	return sqrt(2 / M_PI) * pow(1 / (boltzmannСonstant * temperature), 1.5) * ((-1) * boltzmannСonstant * temperature * v * exp((-1) * v * v / (2 * boltzmannСonstant * temperature)) + pow(boltzmannСonstant * temperature, 1.5) * sqrt(M_PI / 2) * erf(v / sqrt(2 * boltzmannСonstant * temperature)));
}

void CalculateAndRecordTheoreticalMaxwellDistribution() {
	for (int i = 0; i < intervalsNumber; ++i)
		foutTheorMaxwell << speedInterval * i << "\t" << IntegratedMaxwellDistribution(temperature, speedInterval * (i + 1)) - IntegratedMaxwellDistribution(temperature, speedInterval * i) << endl;
}

void RecordNormalizedMaxwellDistribution() {
	int intervalsSum = 0;

	for (int i = 0; i < intervalsNumber; ++i)
		intervalsSum += intervals[i];

	for (int i = 0; i < intervalsNumber; ++i)
		foutMaxwell << double(intervals[i]) / intervalsSum << endl;
}

void main() {
	int currentStep = 0;	// номер текущей итерации

	Initialization();

	for (currentTime = 0; currentTime <= calculationTime; currentTime += timeStep) {
		CalculateEnergyAndForces();
		SolveMotionEquation();
		if (currentStep % 10 == 0)	// запись в файл каждые 10 шагов
			RecordState();
		if (currentStep % 1000 == 0)
			cout << "\t State " << currentStep << " recorded." << endl;
		if (currentTime >= 4  && currentStep % 10 == 0) {
			CalculateMaxwellDistribution();
		}
		++currentStep;
	}

	RecordNormalizedMaxwellDistribution();
	CalculateAndRecordTheoreticalMaxwellDistribution();

	foutEnergy.close();
//	foutPosition.close();
	foutMaxwell.close();
	foutTheorMaxwell.close();

	cout << "Done!" << endl;
}
