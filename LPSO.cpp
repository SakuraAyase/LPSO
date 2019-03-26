#include<cmath>
#include<vector>
#include<ctime>
#include<random>
#include<iostream>
#define pi 3.1415926535
#define E  2.71828182845904523536

using namespace std;

double randDouble(double min, double max)
{
	static default_random_engine engine(time(nullptr));
	uniform_real_distribution<double> dis(min, max);
	return dis(engine);
}
class Particle
{
public:
	double fitness;
	vector<double> position;
	vector<double>velocity;
	vector<double>pBest;
	vector<double>lBest;
	double pBestFitness;
	double lBestFitness;
	Particle() {}

	Particle(vector<double> position, vector<double>velocity, vector<double>best_position, double best_fitness)
	{
		this->position = position;
		this->velocity = velocity;
		this->pBest = best_position;
		this->pBestFitness = best_fitness;
	}
};

bool better(double a, double b)
{
	if (a < b)
		return true;
	else
		return false;
}

class PSO
{
public:

	PSO(int dim, int m, int Tmax, double max, double min, double c1, double c2, double wmax, double wmin, double dt, double percent)
	{
		this->dim = dim;
		this->m = m;
		this->Tmax = Tmax;
		this->max = max;
		this->min = min;
		this->c1 = c1;
		this->c2 = c2;
		this->wmax = wmax;
		this->wmin = wmin;
		this->dt = dt;
		this->percent = percent;
		particles.resize(m);
		
	}


	double fitnessFunction(vector<double> pos)
	{
		double result = 0.0;
		for (int i = 0; i < dim; i++)
		{
			result += pow(pos[i], 2);
		}

		/*
		double temp1 = 0.0;
		double temp2 = 0.0;
		for (int i = 0; i < dim; i++)
		{
			double x = particle.position[i];
			temp1 += pow(x, 2);
			temp2 += cos(2 * pi*x);
		}

		result = -20 * exp(-0.2*sqrt(temp1 / dim)) + 20 + E - exp(temp2 / dim);
		*/
		return result;
	}

	void initialParticles(int i)
	{
		particles[i].position.resize(dim);
		particles[i].velocity.resize(dim);
		particles[i].pBest.resize(dim);
		particles[i].lBest.resize(dim);
		for (int j = 0; j < dim; j++)
		{
			double range = percent * (max - min);
			particles[i].position[j] = randDouble(this->min, this->max);
			particles[i].velocity[j] = randDouble(-range, range);
			particles[i].pBest[j] = particles[i].position[j];
			particles[i].lBest[j] = particles[i].position[j];
		}
		particles[i].fitness = fitnessFunction(particles[i].position);
		particles[i].pBestFitness = fitnessFunction(particles[i].pBest);
		particles[i].lBestFitness = fitnessFunction(particles[i].lBest);
	}

	void initialAllParticles()
	{

		for (int i = 0; i < m; i++)
		{
			initialParticles(i);
		}
	}

	void inertiaWeight()
	{
		//w = randDouble(0.4, 0.6);
		double t = T / ((double)Tmax);
		w = wmax - (wmax - wmin)*t;
	}

	void updateParticle(int i)
	{
		int best = i - 1;
		if (i - 1 < 0)
			best = m - 1;
		if (fitnessFunction(particles[best].pBest) > fitnessFunction(particles[(i + 1) % m].pBest))
		{
			best = (i + 1) % m;
		}
		particles[i].lBest = particles[best].pBest;

		for (int j = 0; j < dim; j++)
		{
			double last_position = particles[i].position[j];
			double range = percent * (max - min);

			particles[i].velocity[j] = w * particles[i].velocity[j] +
				c1 * randDouble(0, 1) * (particles[i].pBest[j] - particles[i].position[j])
				+ c2 * randDouble(0, 1) * (particles[i].lBest[j] - particles[i].position[j]);
			particles[i].position[j] += dt * particles[i].velocity[j];

			if (particles[i].velocity[j] > range)
				particles[i].velocity[j] = range;

			if (particles[i].velocity[j] < -range)
			{
				particles[i].velocity[j] = -range;
			}

			if (particles[i].position[j] > max)
			{
				double thre = randDouble(0, 1);
				if (last_position == max)
				{
					particles[i].position[j] = randDouble(min, max);
				}
				else if (thre < 0.5)
				{
					particles[i].position[j] = max - (max - last_position) * randDouble(0, 1);
				}
				else
				{
					particles[i].position[j] = max;
				}
			}
			if (particles[i].position[j] < min)
			{
				double thre = randDouble(0, 1);
				if (last_position == min)
				{
					particles[i].position[j] = randDouble(min, max);
				}
				else if (thre < 0.5)
				{
					particles[i].position[j] = min + (last_position - min) * randDouble(0, 1);
				}
				else

					particles[i].position[j] = min;
			}

		}
		particles[i].fitness = fitnessFunction(particles[i].position);

		if (particles[i].fitness < particles[i].pBestFitness)
		{
			particles[i].pBestFitness = particles[i].fitness;
			for (int j = 0; j < dim; j++)
			{
				particles[i].pBest[j] = particles[i].position[j];
			}
		}

	}


	void updateAllParticles()
	{
		inertiaWeight();
		for (int i = 0; i < m; i++)
		{
			updateParticle(i);
		}
		T++;
	}

	double getFitness()
	{
		int index = 0;
		for (int i = 0; i < m; i++)
		{
			if (fitnessFunction(particles[i].pBest) < fitnessFunction(particles[index].pBest))
				index = i;
		}
		return fitnessFunction(particles[index].pBest);
	}
private:
	int dim;
	int m;//number of instances

	int T;
	int Tmax;

	double w;
	double max;
	double min;
	double c1;
	double c2;
	double wmax;
	double wmin;

	double dt;//时间步长
	double percent;


	vector<Particle> particles;


};

void run(vector<double>& result1)
{
	int dim = 30;
	int m = 20;
	int Tmax = 2000;
	double max = 100;
	double min = -100;
	double c1 = 2;
	double c2 = 2;
	double wmax = 0.9;
	double wmin = 0.4;
	double dt = 1.0;
	double percent = 0.2;

	PSO pso = PSO(dim, m, Tmax, max, min, c1, c2, wmax, wmin, dt, percent);
	pso.initialAllParticles();

	vector<double>fitness;
	fitness.push_back(pso.getFitness());

	for (int i = 0; i < Tmax; i++)
	{
		pso.updateAllParticles();
		cout << ":";
		//fitness.push_back(pso.getFitness());
		fitness.push_back(pso.getFitness());
		cout << "第" << i << "次迭代结果：";
		cout << ", fitness = " << pso.getFitness() << endl;
	}

	result1 = fitness;
}

int main()
{

	int times = 5;
	int interval = 10;
	vector<double> result1;

	run(result1);

	for (int i = 1; i < times; i++)
	{
		vector<double> result1_temp;
		run(result1_temp);
		for (int j = 0; j < result1_temp.size(); j++)
		{
			result1[j] += result1_temp[j];
		}
	}
	for (int j = 0; j < result1.size(); j++)
	{
		result1[j] /= times;
	}

	for (int j = 0; j < result1.size(); j++)
	{
		if (j%interval == 0)
			cout << result1[j] << " ";
	}

	system("pause");
}