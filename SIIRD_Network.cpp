#include <iostream>
#include <random>
#include <cstdlib>
#include <math.h>
#include <vector>
#include <map>
#include <chrono>
#include <string>
#include <fstream>
#include <filesystem>

// A matrix struct for convenient head memory allocation
template <typename T>
struct Matrix
{
	T** dptr = nullptr; // Double pointer which is the main part of handling values
	T* ptr = nullptr; // Regular pointer to ensure that memory is allocated in one chunk
	int rows = 0; // Amount of rows in the matrix
	int columns = 0; // Amount of columns in the matrix

	// Returns the head allocated double pointer
	T** getDptr()
	{
		return dptr;
	}

	// Set the dimensions of matrix
	void setDim(int A, int B)
	{
		if (dptr || ptr)
		{
			delete[] dptr;
			delete[] ptr;
		}

		dptr = (T**)malloc(A * sizeof(T*));
		ptr = (T*)malloc(A * B * sizeof(T));
		rows = A;
		columns = B;

		for (int i = 0; i < A; i++)
			dptr[i] = ptr + i * B;
	}

	// Initialize all values of matrix to 0
	void matrixInit()
	{
		for (int x = 0; x < rows; x++)
		{
			for (int y = 0; y < columns; y++)
			{
				dptr[x][y] = 0.0;
			}
		}
	}

	// Set all diagonal values to 0
	void zerosDiagonal()
	{
		if (rows == columns)
		{
			for (int i = 0; i < rows; i++)
				dptr[i][i] = 0;
		}
		else
			std::cout << "Matrix needs to have same amount of rows and columns\n";
	}

	// Prints the dimensions of matrix
	void getDim()
	{
		if (rows || columns)
			std::cout << "Rows: " << rows << ", Columns: " << columns << "\n";
		else
			std::cout << "Dimensions haven't been set\n";
	}

	// Displays a spesific element of the matrix
	void displayElement()
	{
		std::cout << "Dimensions of matrix are " << rows << " x " << columns << "\n";
		std::cout << "Which element would you like to check?\n";
		std::cout << "Give a negative value to stop\n";
		while (true)
		{
			int x;
			int y;
			std::cin >> x >> y;
			if (x < 0 || y < 0)
				break;
			if (x > rows - 1 || y > columns - 1)
				std::cout << "Values are too big\n";
			else
				std::cout << "Value found " << (float)dptr[x][y] << "\n";
		}
	}

	// Displays a spesific submatrix of the original matrix
	void displaySubmatrix()
	{
		std::cout << "Dimensions of matrix are " << rows << " x " << columns << "\n";
		std::cout << "Which submatrix would you like to check?\n";
		std::cout << "Give ranges in form: x_1 x_2 y_1 y_2\n";
		std::cout << "Give a negative value to stop\n";

		while (true)
		{
			int x_1;
			int x_2;
			int y_1;
			int y_2;
			std::cin >> x_1 >> x_2 >> y_1 >> y_2;
			if (x_1 < 0 || x_2 < 0 || y_1 < 0 || y_2 < 0)
				break;
			if (x_1 > x_2 || y_1 > y_2)
				std::cout << "Incorrect range\n";
			else if (x_1 > rows - 1 || x_2 > rows - 1 || y_1 > columns - 1 || y_2 > columns - 1)
				std::cout << "Values are out of range\n";
			else
			{
				for (int i = x_1; i <= x_2; i++)
				{
					for (int j = y_1; j <= y_2; j++)
					{
						std::cout << dptr[i][j] << " ";
					}
					std::cout << "\n";
				}	
				std::cout << "\n";
			}
		}
	}

	Matrix()
	{
	}

	Matrix(int A, int B)
	{
		dptr = (T**)malloc(A * sizeof(T*));
		ptr = (T*)malloc(A * B * sizeof(T));
		rows = A;
		columns = B;

		for (int i = 0; i < A; i++)
			dptr[i] = ptr + i * B;
	}

	~Matrix()
	{
		if (dptr || ptr)
		{
			delete[] dptr;
			delete[] ptr;
		}
	}
};

// Network class for simulating SIIRD epidemic model
class Network_SIIRD
{
private:
	int GridSize = 10000; // Dimensions of square where coordinates can generate
	int nodes_ = 0; 
	int* Nodecount = &nodes_; // Amount of nodes in the network

	unsigned long generator_seed;
	std::default_random_engine generator_;
	std::default_random_engine* Generator = &generator_; // Generator for pseudo random values

	Matrix<float> adjacency_;
	float** AdjacencyMatrix; // Adjacencymatrix where all information about network is stored

	Matrix<float> familyadjacency_;
	float** FamilyAdjacencyMatrix; // Adjacencymatrix where information about families in network is stored

	bool* Susceptiple = nullptr; // Pointer which stores information about susceptible individuals
	int* Infected1 = nullptr; // Pointer which stores information about individuals in infected group 1
	int* bufferinf1_ = nullptr; // Buffer pointer for infected group 1
	int* Infected2 = nullptr; // Pointer which stores information about individuals in infected group 2
	int* bufferinf2_ = nullptr; // Buffer pointer for infected group 2
	bool* Recovered = nullptr; // Pointer which stores information about recovered individuals
	bool* Deceased = nullptr; // Pointer which stores information about deceased individuals

	int suscount_;
	int* AmountOfSusceptiple = & suscount_; // Current amount of susceptible individuals
	int inf1count_;
	int* AmountOfInfected1 = &inf1count_; // Current amount of individuals in infected group 1
	int inf2count_;
	int* AmountOfInfected2 = &inf2count_; // Current amount of individuals in infected group 2
	int reccount_;
	int* AmountOfRecovered = &reccount_; // Current amount of recovered individuals
	int deccount_;
	int* AmountOfDeceased = &deccount_; // Current amount of deceased individuals

	const float DeathRate = 0.01; // Death rate of epidemy

public:
	enum Age { adult = 1, upper = 2, pre = 3 };

private:

	// Random value between two values from a uniform distribution
	int RandomUniform(int first_number, int last_number, std::default_random_engine* generator)
	{
		std::uniform_int_distribution<> uniform_distribution(first_number, last_number);
		int random_number = uniform_distribution(*generator);
		return random_number;
	}

	// Random value from a binomial distribution
	int RandomBinom(int how_many, double probability, std::default_random_engine* generator)
	{
		std::binomial_distribution<int> binom_distribution(how_many, probability);
		int random_number = binom_distribution(*generator);
		return random_number;
	}

	// Calculates the shortest euclidean distance to household
	int ShortestDistance(int house_x, int house_y, int amount_of_mixing_places, int** mix_coordinates)
	{
		double distance1 = NULL;
		int chosenmix_ = NULL;

		// Check the location of all mixing places and choose the one which is the closest to the house
		for (int i = 0; i < amount_of_mixing_places; i++)
		{
			double distance2 = sqrt(pow(house_x - mix_coordinates[i][0], 2) + pow(house_y - mix_coordinates[i][1], 2));
			if (distance1 == NULL || distance2 < distance1)
			{
				distance1 = distance2;
				chosenmix_ = i;
			}
		}
		return chosenmix_;
	}

	// Updates the coefficient in an adjacency matrix
	float UpdateCoefficient(float old_value, float new_value)
	{
		float updated_ = 1 - ((1 - old_value) * (1 - new_value));
		return updated_;
	}

	// Updates the adjacency matrix
	void UpdateAdjacency(float new_value, std::map<int, std::vector<int>> discriminating_values, float** adjacency_matrix)
	{
		// Loop through different mixing places and the values assigned to each one
		for (const std::pair<const int, std::vector<int>>& values : discriminating_values)
		{
			for (int value_x : values.second)
			{
				for (int value_y : values.second)
				{
					if (value_x != value_y)
						adjacency_matrix[value_x][value_y] = UpdateCoefficient(adjacency_matrix[value_x][value_y], new_value);
				}
			}
		}
	}

	// Function which determines if an indected individual stays home
	bool StayHome(int& inf1, int& inf2)
	{
		float stayhome_ = ((float)(8 - inf1)) / 10.0;
		if (stayhome_ > ((float)rand()) / RAND_MAX)
		{
			inf2 = inf1;
			inf1 = 0;
			return 1;
		}
		return 0;
	}

	// Initializes pointers
	void nodeInit()
	{
		std::uniform_int_distribution<> uniform_distribution(0, *Nodecount);

		// Initialize each node a susceptible
		for (int i = 0; i < *Nodecount; i++)
		{
			Susceptiple[i] = 1;
			Infected1[i] = 0;
			Infected2[i] = 0;
			bufferinf1_[i] = 0;
			bufferinf2_[i] = 0;
			Recovered[i] = 0;
			Deceased[i] = 0;
		}

		// Assign initial infected group 1 individuals
		for (int i = 0; i < *AmountOfInfected1; i++)
		{
			while (true)
			{
				int node_ = uniform_distribution(*Generator);
				if (Susceptiple[node_] == 1)
				{
					Infected1[node_] = 7;
					Susceptiple[node_] = 0;
					break;
				}
			}
		}

		// Assign initial infected group 2 individuals
		for (int i = 0; i < *AmountOfInfected2; i++)
		{
			while (true)
			{
				int node_ = uniform_distribution(*Generator);
				if (Susceptiple[node_] == 1)
				{
					Infected2[node_] = 7;
					Susceptiple[node_] = 0;
					break;
				}
			}
		}

		// Assign initial recovered individuals
		for (int i = 0; i < *AmountOfRecovered; i++)
		{
			while (true)
			{
				int node_ = uniform_distribution(*Generator);
				if (Susceptiple[node_] == 1)
				{
					Recovered[node_] = 7;
					Susceptiple[node_] = 0;
					break;
				}
			}
		}
	}

	// Function to infect random people
	void InfectRandom()
	{
		int randomnode_ = RandomUniform(0, *Nodecount, Generator);
		if (0.001 > ((double)rand()) / RAND_MAX && Susceptiple[randomnode_])
		{
			bufferinf1_[randomnode_] = 7;
			Susceptiple[randomnode_] = 0;

			*AmountOfInfected1 += 1;
			*AmountOfSusceptiple -= 1;
		}
	}

	// Function that simulates the next days happenings
	void InfectSusceptiple()
	{
		for (int i = 0; i < *Nodecount; i++)
		{
			// Check if infected is in infected group 1
			if (Infected1[i])
			{
				for (int j = 0; j < *Nodecount; j++)
				{
					if (AdjacencyMatrix[i][j] > ((double)rand()) / RAND_MAX && Susceptiple[j])
					{
						bufferinf1_[j] = 7;
						Susceptiple[j] = 0;

						*AmountOfInfected1 += 1;
						*AmountOfSusceptiple -= 1;
					}
				}

				InfectRandom();

				if (DeathRate > ((double)rand()) / RAND_MAX)
				{
					Deceased[i] = 1;
					Infected1[i] = 0;
					bufferinf1_[i] = 0;

					*AmountOfInfected1 -= 1;
					*AmountOfDeceased += 1;
				}
				else if (StayHome(Infected1[i], bufferinf2_[i]))
				{
					bufferinf1_[i] = 0;
					*AmountOfInfected1 -= 1;
					*AmountOfInfected2 += 1;
				}
				else if ((Infected1[i] -= 1) == 0)
				{
					Recovered[i] = 1;
					bufferinf1_[i] = 0;

					*AmountOfInfected1 -= 1;
					*AmountOfRecovered += 1;
				}
				else
				{
					bufferinf1_[i] = Infected1[i];
				}
			}
			// Check if infected is in infected group 1
			else if (Infected2[i])
			{
				for (int j = 0; j < *Nodecount; j++)
				{
					if (FamilyAdjacencyMatrix[i][j] > ((double)rand()) / RAND_MAX && Susceptiple[j])
					{
						bufferinf1_[j] = 7;
						Susceptiple[j] = 0;

						*AmountOfInfected1 += 1;
						*AmountOfSusceptiple -= 1;
					}
				}

				if (DeathRate > ((double)rand()) / RAND_MAX)
				{
					Deceased[i] = 1;
					Infected2[i] = 0;
					bufferinf2_[i] = 0;

					*AmountOfInfected2 -= 1;
					*AmountOfDeceased += 1;
				}
				else if ((Infected2[i] -= 1) == 0)
				{
					Recovered[i] = 1;
					bufferinf2_[i] = 0;

					*AmountOfInfected2 -= 1;
					*AmountOfRecovered += 1;
				}
				else
				{
					bufferinf2_[i] = Infected2[i];
				}
			}
		}

		memcpy(Infected1, bufferinf1_, (*Nodecount) * sizeof(int));
		memcpy(Infected2, bufferinf2_, (*Nodecount) * sizeof(int));
	}

public:

	// Generates households and assings each house a family of size 1+
	void GenerateFamilies(float mixing_rate, int amount_of_houses, int** house_configuration, int** house_coordinates)
	{
		std::normal_distribution<double> normal_distribution(2.8, 1.5);

		for (int i = 0; i < amount_of_houses; i++)
		{
			int family_size = round(normal_distribution(*Generator));
			if (family_size <= 0)
				*Nodecount += 1;
			else
				*Nodecount += family_size;
			house_configuration[i][0] = family_size;

			if (family_size <= 1)
				house_configuration[i][1] = 1, house_configuration[i][2] = 0, house_configuration[i][3] = 0;
			else if (family_size == 2)
				house_configuration[i][1] = 2, house_configuration[i][2] = 0, house_configuration[i][3] = 0;
			else
			{
				int old_children = RandomBinom(family_size - 2, 2.0/3.0, Generator);
				int young_children = family_size - 2 - old_children;
				house_configuration[i][1] = 2, house_configuration[i][2] = old_children, house_configuration[i][3] = young_children;
			}

			house_coordinates[i][0] = RandomUniform(0, 10000, Generator), house_coordinates[i][1] = RandomUniform(0, 10000, Generator);

		}

		adjacency_.setDim(*Nodecount, *Nodecount);
		familyadjacency_.setDim(*Nodecount, *Nodecount);
		adjacency_.matrixInit();
		familyadjacency_.matrixInit();
		AdjacencyMatrix = adjacency_.getDptr();
		FamilyAdjacencyMatrix = familyadjacency_.getDptr();

		// Initialize household connections to adjacencymatrices
		int nextnode_ = 0;
		for (int household = 0; household < amount_of_houses; household++)
		{
			int current_size = house_configuration[household][0];
			for (int node_x = 0; node_x < current_size; node_x++)
			{
				for (int node_y = 0; node_y < current_size; node_y++)
				{
					if (node_x != node_y)
					{
						AdjacencyMatrix[nextnode_ + node_x][nextnode_ + node_y] = mixing_rate;
						FamilyAdjacencyMatrix[nextnode_ + node_x][nextnode_ + node_y] = mixing_rate;
					}
				}
			}

			nextnode_ += current_size;
		}
	}

	// Generates mixing places for disease and assings each household randomly to them
	void MixingRandom(float mixing_rate, int amount_of_mixing_places, int amount_of_homes, int** house_configuration, Network_SIIRD::Age age)
	{
		std::map<int, std::vector<int>> discriminator_;
		std::uniform_int_distribution<> uniform_distribution(0, amount_of_mixing_places - 1);

		// Grouping elements
		int nextnode_ = 0;
		if (age == 1)
		{
			for (int i = 0; i < amount_of_homes; ++i)
			{
				int randomvalue_ = uniform_distribution(*Generator);
				for (int j = 0; j < house_configuration[i][1]; j++)
					discriminator_[randomvalue_].push_back(nextnode_ + j);

				nextnode_ += house_configuration[i][0];
			}
		}
		else if (age == 2)
		{
			for (int i = 0; i < amount_of_homes; ++i)
			{
				nextnode_ += house_configuration[i][1];
				int randomvalue_ = uniform_distribution(*Generator);

				for (int j = 0; j < house_configuration[i][2]; j++)
				{
					discriminator_[randomvalue_].push_back(nextnode_);
					nextnode_++;
				}

				nextnode_ += house_configuration[i][3];
			}
		}
		else if (age == 3)
		{
			for (int i = 0; i < amount_of_homes; ++i)
			{
				nextnode_ = nextnode_ + house_configuration[i][1] + house_configuration[i][2];
				int randomvalue_ = uniform_distribution(*Generator);

				for (int j = 0; j < house_configuration[i][3]; j++)
				{
					discriminator_[randomvalue_].push_back(nextnode_);
					nextnode_++;
				}
			}
		}
		
		UpdateAdjacency(mixing_rate, discriminator_, AdjacencyMatrix);
	}

	// Generates mixing places for disease and assings each household based on shortest distance to them
	void MixingDistance(float mixing_rate, int amount_of_mixing_places, int amount_of_homes, int** house_configuration, int** house_coordinates, Network_SIIRD::Age age)
	{
		std::map<int, std::vector<int>> discriminator_;

		Matrix<int> m_mix_coordinates(amount_of_mixing_places, 2);
		int** mixing_coordinates = m_mix_coordinates.getDptr();

		for (int i = 0; i < amount_of_mixing_places; i++)
			mixing_coordinates[i][0] = RandomUniform(0, 10000, Generator), mixing_coordinates[i][1] = RandomUniform(0, 10000, Generator);

		// Grouping elements
		int nextnode_ = 0;
		if (age == 1)
		{
			for (int i = 0; i < amount_of_homes; ++i)
			{
				int chosenvalue_ = ShortestDistance(house_coordinates[i][0], house_coordinates[i][1], amount_of_mixing_places, mixing_coordinates);

				for (int j = 0; j < house_configuration[i][1]; j++)
					discriminator_[chosenvalue_].push_back(nextnode_ + j);

				nextnode_ += house_configuration[i][0];
			}
		}
		else if (age == 2)
		{
			for (int i = 0; i < amount_of_homes; ++i)
			{
				nextnode_ += house_configuration[i][1];
				int chosenvalue_ = ShortestDistance(house_coordinates[i][0], house_coordinates[i][1], amount_of_mixing_places, mixing_coordinates);

				for (int j = 0; j < house_configuration[i][2]; j++)
				{
					discriminator_[chosenvalue_].push_back(nextnode_);
					nextnode_++;
				}

				nextnode_ += house_configuration[i][3];
			}
		}
		else if (age == 3)
		{
			for (int i = 0; i < amount_of_homes; ++i)
			{
				nextnode_ = nextnode_ + house_configuration[i][1] + house_configuration[i][2];
				int chosenvalue_ = ShortestDistance(house_coordinates[i][0], house_coordinates[i][1], amount_of_mixing_places, mixing_coordinates);

				for (int j = 0; j < house_configuration[i][3]; j++)
				{
					discriminator_[chosenvalue_].push_back(nextnode_);
					nextnode_++;
				}
			}
		}

		UpdateAdjacency(mixing_rate, discriminator_, AdjacencyMatrix);
	}

	// Allocate head memory for pointers
	void allocPtrs()
	{
		Susceptiple = (bool*)malloc(*Nodecount * sizeof(bool));
		Infected1 = (int*)malloc(*Nodecount * sizeof(int));
		bufferinf1_ = (int*)malloc(*Nodecount * sizeof(int));
		Infected2 = (int*)malloc(*Nodecount * sizeof(int));
		bufferinf2_ = (int*)malloc(*Nodecount * sizeof(int));
		Recovered = (bool*)malloc(*Nodecount * sizeof(bool));
		Deceased = (bool*)malloc(*Nodecount * sizeof(bool));
	}

	// Function which simulates days of the disease until no infected people are left

	void SimulateDisease(int numberofsimulation_)
	{
		suscount_ = 0;
		inf1count_ = 0;
		inf2count_ = 0;
		reccount_ = 0;
		deccount_ = 0;

		int daysofepidemy_ = 0;

		std::string filename_ = "Simulation " + std::to_string(numberofsimulation_) + ".txt";

		std::ofstream Diseasedata(filename_);

		std::cout << "There are " << *Nodecount << " nodes in the network.\n";
		std::cout << "How many initial infectees, i.e. infected1 and infected2 and recovered?\n";
		std::cout << "Give answer in form: inf1 inf2 rec\n";

		while (true)
		{
			std::cin >> *AmountOfInfected1 >> *AmountOfInfected2 >> *AmountOfRecovered;
			if (*AmountOfInfected1 + *AmountOfInfected2 + *AmountOfRecovered > *Nodecount)
			{
				std::cout << "Too many, can't be greater than " << *Nodecount << "\n";
			}
			else
				break;
		}

		*AmountOfSusceptiple = *Nodecount - *AmountOfInfected1 - *AmountOfInfected2 - *AmountOfRecovered;

		nodeInit();
		Diseasedata << *AmountOfSusceptiple << ";" << *AmountOfInfected1 << ";" << *AmountOfInfected2 << ";" << *AmountOfRecovered << ";" << *AmountOfDeceased << "\n";

		while (inf1count_ || inf2count_)
		{
			std::cout << "This is day " << daysofepidemy_ << " of epidemic\n";
			daysofepidemy_ += 1;
			InfectSusceptiple();
			Diseasedata << *AmountOfSusceptiple << ";" << *AmountOfInfected1 << ";" << *AmountOfInfected2 << ";" << *AmountOfRecovered << ";" << *AmountOfDeceased << "\n";
		}

		Diseasedata << "Days of epidemic: " << daysofepidemy_;
		Diseasedata.close();
	}

	// Variation of the same function, which is adjusted for continous simulation
	void SimulateDisease(int numberofsimulation_, int householdsize_, int inf1init_, int inf2init_, float recoveredproportion_, std::string filepath_)
	{
		suscount_ = 0;
		inf1count_ = 0;
		inf2count_ = 0;
		reccount_ = 0;
		deccount_ = 0;

		std::string filename_ = "Simulation_" + std::to_string(numberofsimulation_) + "_" + std::to_string(householdsize_) + "_" + std::to_string(inf1init_) + "_" + std::to_string(inf2init_) + "_" + std::to_string(recoveredproportion_) + ".txt";

		std::ofstream Diseasedata(filepath_ + filename_);

		*AmountOfInfected1 = inf1init_;
		*AmountOfInfected2 = inf2init_;
		*AmountOfRecovered = *Nodecount * recoveredproportion_;
		*AmountOfSusceptiple = *Nodecount - *AmountOfInfected1 - *AmountOfInfected2 - *AmountOfRecovered;

		nodeInit();
		Diseasedata << (float)*AmountOfSusceptiple / (float)*Nodecount << ";" << (float)*AmountOfInfected1 / (float)*Nodecount << ";" << (float)*AmountOfInfected2 / (float)*Nodecount << ";" << (float)*AmountOfRecovered / (float)*Nodecount << ";" << (float)*AmountOfDeceased / (float)*Nodecount << "\n";

		while (inf1count_ || inf2count_)
		{
			std::cout << "This is day " << daysofepidemy_ << " of epidemic\n";
			InfectSusceptiple();
			Diseasedata << (float)*AmountOfSusceptiple / (float)*Nodecount << ";" << (float)*AmountOfInfected1 / (float)*Nodecount << ";" << (float)*AmountOfInfected2 / (float)*Nodecount << ";" << (float)*AmountOfRecovered / (float)*Nodecount << ";" << (float)*AmountOfDeceased / (float)*Nodecount << "\n";
		}

		Diseasedata.close();
	}

	// Function which prints out the nodecount
	void getNetworksize()
	{
		std::cout << "Network has " << *Nodecount << " nodes\n";
	}

	// Displays a spesific element of the adjacencymatrix
	void displayElement()
	{
		adjacency_.displayElement();
	}

	// Displays a spesific submatrix of the adjacencymatrix
	void displaySubmatrix()
	{
		adjacency_.displaySubmatrix();
	}

	Network_SIIRD(unsigned long seed)
	{
		generator_seed = seed;
		(*Generator).seed(generator_seed);
	}

	~Network_SIIRD()
	{
		if (Susceptiple || Infected1 || bufferinf1_ || Infected2 || bufferinf2_ || Recovered || Deceased)
		{
			free(Susceptiple);
			free(Infected1);
			free(bufferinf1_);
			free(Infected2);
			free(bufferinf2_);
			free(Recovered);
			free(Deceased);
		}
	}
};

int main()
{
	int house_;
	int* n_House = &house_;
	int mixing_Pre;
	int* n_Pre = &mixing_Pre;
	int mixing_Upper;
	int* n_Upper = &mixing_Upper;
	int mixing_Work;
	int* n_Work = &mixing_Work;

	std::cout << "How many households, preschools, upper schools and workplaces? Separate each with \" \":\n";
	std::cin >> house_ >> mixing_Pre >> mixing_Upper >> mixing_Work;

	Matrix<int> M_house_conf(*n_House, 4);
	Matrix<int> M_house_coor(*n_House, 2);
	int** House_Conf = M_house_conf.getDptr();
	int** House_Coor = M_house_coor.getDptr();

	auto start = std::chrono::high_resolution_clock::now();

	Network_SIIRD Disease(1);
	Disease.GenerateFamilies(0.005, *n_House, House_Conf, House_Coor);
	Disease.MixingRandom(0.001, *n_Work, *n_House, House_Conf, Network_SIIRD::adult);
	Disease.MixingDistance(0.003, *n_Upper, *n_House, House_Conf, House_Coor, Network_SIIRD::upper);
	Disease.MixingDistance(0.003, *n_Pre, *n_House, House_Conf, House_Coor, Network_SIIRD::pre);
	Disease.allocPtrs();
	Disease.SimulateDisease(0);

	auto end = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
	std::cout << "Simulation took " << duration.count() << " microseconds.\n";
}


