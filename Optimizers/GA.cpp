//****************************************************************************
//
//    A Simple Genetic Algorithm
//
//    This is a simple genetic algorithm implementation where the 
//    evaluation function takes positive values only and the 
//    fitness of an individual is the same as the value of the 
//    objective function.  
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 April 2014
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//	  Modified version by Sohail R. Reddy
//
//  How to Use:
//
//    Initialize the GA using the function 
//
//		void initialize ( int npop, vector<double> lbound, vector<double> ubound, 
//					 double (*fitnessFunc)(vector<double> x), int &tmp_seed );
//
//		Inputs:
//			npop: is the number of members in population
//			lbound: is a vector of size nVar with lower bounds of variables
//			ubound: is a vector of size nVar with lower bounds of variables
//			*fitnessFunc: is a pointer to a function used to compute the fitness given vector 'x'
//
//
//	 Run GA using the function 
//
//		vector<double> runGA (int maxgens, double pmutation, double pxover );
//
//		Inputs:
//			maxgens: is the maximum number of generations
//			pmutation: is the probability of mutation
//			pxover: is the probability of crossover
//	
//		Outputs:
//			vector of size nVar contains the optimum solution
//
//****************************************************************************


# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>
# include <vector>
using namespace std;



int POPSIZE = 0;
int NVARS = 0;
int seed = 0;

struct genotype {
  vector<double> gene;
  vector<double> upper;
  vector<double> lower;
  double fitness;
  double rfitness;
  double cfitness;
}; 
vector<genotype> population;
vector<genotype> newpopulation;



vector<double> runGA (int maxgens, double pmutation, double pxover );
void crossover ( int &seed, double pxover);
void elitist ( );
void evaluate ();
int i4_uniform_ab ( int a, int b, int &seed );
void initialize ( int npop, vector<double> lbound, vector<double> ubound, double (*fitnessFunc)(vector<double> x), int &in_seed );
void keep_the_best ( );
void mutate ( int &seed, double pmutation );
double r8_uniform_ab ( double a, double b, int &seed );
void report ( int generation );
void selector ( int &seed );
void Xover ( int one, int two, int &seed );
double (*evaluateMember)(vector<double> x);


//****************************************************************************80

vector<double> runGA (int maxgens, double pmutation, double pxover){
//    maxgens is the maximum number of generations.
//    pmutation is the probability of mutation.
//    pxover is the probability of crossover.                          

  int generation;
  int i;
 

  evaluate ( );
  keep_the_best ( );

  for ( generation = 0; generation < maxgens; generation++ )
  {
    selector ( seed );
    crossover ( seed,pxover );
    mutate ( seed , pmutation);
    report ( generation );
    evaluate ( );
    elitist ( );
  }
  
  vector<double> best = population[POPSIZE].gene;
  population.clear();
  newpopulation.clear();
  return best;
}
//****************************************************************************80


void evaluate ()
{
	int member;
	for ( member = 0; member < POPSIZE; member++ )
	{
		population[member].fitness = evaluateMember(population[member].gene);
	 }
	return;
}
//****************************************************************************80


void crossover ( int &seed, double pxover ) {
  const double a = 0.0;
  const double b = 1.0;
  int mem;
  int one;
  int first = 0;
  double x;

  for ( mem = 0; mem < POPSIZE; ++mem )
  {
    x = r8_uniform_ab ( a, b, seed );

    if ( x < pxover )
    {
      ++first;

      if ( first % 2 == 0 )
      {
        Xover ( one, mem, seed );
      }
      else
      {
        one = mem;
      }

    }
  }
  return;
}

//****************************************************************************80

void elitist ( ) {
  int i;
  double best;
  int best_mem;
  double worst;
  int worst_mem;

  best = population[0].fitness;
  worst = population[0].fitness;

  for ( i = 0; i < POPSIZE - 1; ++i )
  {
    if ( population[i+1].fitness < population[i].fitness )
    {

      if ( best <= population[i].fitness )
      {
        best = population[i].fitness;
        best_mem = i;
      }

      if ( population[i+1].fitness <= worst )
      {
        worst = population[i+1].fitness;
        worst_mem = i + 1;
      }

    }
    else
    {

      if ( population[i].fitness <= worst )
      {
        worst = population[i].fitness;
        worst_mem = i;
      }

      if ( best <= population[i+1].fitness )
      {
        best = population[i+1].fitness;
        best_mem = i + 1;
      }

    }

  }
// 
//  If the best individual from the new population is better than 
//  the best individual from the previous population, then 
//  copy the best from the new population; else replace the 
//  worst individual from the current population with the 
//  best one from the previous generation                     
//
  if ( population[POPSIZE].fitness <= best )
  {
    for ( i = 0; i < NVARS; i++ )
    {
      population[POPSIZE].gene[i] = population[best_mem].gene[i];
    }
    population[POPSIZE].fitness = population[best_mem].fitness;
  }
  else
  {
    for ( i = 0; i < NVARS; i++ )
    {
      population[worst_mem].gene[i] = population[POPSIZE].gene[i];
    }
    population[worst_mem].fitness = population[POPSIZE].fitness;
  } 

  return;
}
//****************************************************************************80


int i4_uniform_ab ( int a, int b, int &seed ) {
  int c;
  const int i4_huge = 2147483647;
  int k;
  float r;
  int value;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "I4_UNIFORM_AB - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }
//
//  Guarantee A <= B.
//
  if ( b < a )
  {
    c = a;
    a = b;
    b = c;
  }

  k = seed / 127773;

  seed = 16807 * ( seed - k * 127773 ) - k * 2836;

  if ( seed < 0 )
  {
    seed = seed + i4_huge;
  }

  r = ( float ) ( seed ) * 4.656612875E-10;
//
//  Scale R to lie between A-0.5 and B+0.5.
//
  r = ( 1.0 - r ) * ( ( float ) a - 0.5 ) 
    +         r   * ( ( float ) b + 0.5 );
//
//  Use rounding to convert R to an integer between A and B.
//
  value = round ( r );
//
//  Guarantee A <= VALUE <= B.
//
  if ( value < a )
  {
    value = a;
  }
  if ( b < value )
  {
    value = b;
  }

  return value;
}
//****************************************************************************80

void initialize ( int npop, vector<double> lbound, vector<double> ubound, double (*fitnessFunc)(vector<double> x), int &in_seed ) {
//    npop is the number of members in populaiton.
//    lbound is the lower bound on each variables
//    ubound is the upper bound on each variables
//	  fitnessFunc is the function pointer to compute the fitness of a single member
//	  in_seed is the random seed

  int i;
  int j;

	evaluateMember = fitnessFunc;

	seed = in_seed;
	NVARS = lbound.size();
	POPSIZE = npop;
	
	population.resize(POPSIZE+1);
	newpopulation.resize(POPSIZE+1);
  
    for ( j = 0; j < POPSIZE+1; j++ ) {
		population[j].lower.resize(NVARS);
		population[j].upper.resize(NVARS);
		population[j].gene.resize(NVARS);
		
		newpopulation[j].lower.resize(NVARS);
		newpopulation[j].upper.resize(NVARS);
		newpopulation[j].gene.resize(NVARS);
	}  
  
// 
//  Initialize variables within the bounds 
//
  for ( i = 0; i < NVARS; i++ )
  {
    for ( j = 0; j < POPSIZE; j++ )
    {
      population[j].fitness = 0;
      population[j].rfitness = 0;
      population[j].cfitness = 0;
      population[j].lower[i] = lbound[i];
      population[j].upper[i]= ubound[i];
      population[j].gene[i] = r8_uniform_ab ( lbound[i], ubound[i], seed );
    }
  }

  return;
}




void keep_the_best ( ) {
  int cur_best;
  int mem;
  int i;

  cur_best = 0;

  for ( mem = 0; mem < POPSIZE; mem++ )
  {
    if ( population[POPSIZE].fitness < population[mem].fitness )
    {
      cur_best = mem;
      population[POPSIZE].fitness = population[mem].fitness;
    }
  }
// 
//  Once the best member in the population is found, copy the genes.
//


  for ( i = 0; i < NVARS; i++ )
  {
    population[POPSIZE].gene[i] = population[cur_best].gene[i];
  }

  return;
}
//****************************************************************************80

void mutate ( int &seed, double pmutation ) {
  const double a = 0.0;
  const double b = 1.0;
  int i;
  int j;
  double lbound;
  double ubound;
  double x;

  for ( i = 0; i < POPSIZE; i++ )
  {
    for ( j = 0; j < NVARS; j++ )
    {
      x = r8_uniform_ab ( a, b, seed );
      if ( x < pmutation )
      {
        lbound = population[i].lower[j];
        ubound = population[i].upper[j];  
        population[i].gene[j] = r8_uniform_ab ( lbound, ubound, seed );
      }
    }
  }

  return;
}
//****************************************************************************80

double r8_uniform_ab ( double a, double b, int &seed ) {
  int i4_huge = 2147483647;
  int k;
  double value;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "R8_UNIFORM_AB - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  k = seed / 127773;

  seed = 16807 * ( seed - k * 127773 ) - k * 2836;

  if ( seed < 0 )
  {
    seed = seed + i4_huge;
  }

  value = ( double ) ( seed ) * 4.656612875E-10;

  value = a + ( b - a ) * value;

  return value;
}
//****************************************************************************80

void report ( int generation ) {
  double avg;
  double best_val;
  int i;
  double square_sum;
  double stddev;
  double sum;
  double sum_square;

  if ( generation == 0 )
  {
    cout << "\n";
    cout << "  Generation       Best            Average       Standard \n";
    cout << "  number           value           fitness       deviation \n";
    cout << "\n";
  }

  sum = 0.0;
  sum_square = 0.0;

  for ( i = 0; i < POPSIZE; i++ )
  {
    sum = sum + population[i].fitness;
    sum_square = sum_square + population[i].fitness * population[i].fitness;
  }

  avg = sum / ( double ) POPSIZE;
  square_sum = avg * avg * POPSIZE;
  stddev = sqrt ( ( sum_square - square_sum ) / ( POPSIZE - 1 ) );
  best_val = population[POPSIZE].fitness;

  cout << "  " << setw(8) << generation 
       << "  " << setw(14) << best_val 
       << "  " << setw(14) << avg 
       << "  " << setw(14) << stddev << "\n";

  return;
}
//****************************************************************************80

void selector ( int &seed ) {
  const double a = 0.0;
  const double b = 1.0;
  int i;
  int j;
  int mem;
  double p;
  double sum;
//
//  Find the total fitness of the population.
//
  sum = 0.0;
  for ( mem = 0; mem < POPSIZE; mem++ )
  {
    sum = sum + population[mem].fitness;
  }
//
//  Calculate the relative fitness of each member.
//
  for ( mem = 0; mem < POPSIZE; mem++ )
  {
    population[mem].rfitness = population[mem].fitness / sum;
  }
// 
//  Calculate the cumulative fitness.
//
  population[0].cfitness = population[0].rfitness;
  for ( mem = 1; mem < POPSIZE; mem++ )
  {
    population[mem].cfitness = population[mem-1].cfitness +       
      population[mem].rfitness;
  }
// 
//  Select survivors using cumulative fitness. 
//
  for ( i = 0; i < POPSIZE; i++ )
  { 
    p = r8_uniform_ab ( a, b, seed );
    if ( p < population[0].cfitness )
    {
      newpopulation[i] = population[0];      
    }
    else
    {
      for ( j = 0; j < POPSIZE; j++ )
      { 
        if ( population[j].cfitness <= p && p < population[j+1].cfitness )
        {
          newpopulation[i] = population[j+1];
        }
      }
    }
  }
// 
//  Overwrite the old population with the new one.
//
  for ( i = 0; i < POPSIZE; i++ )
  {
    population[i] = newpopulation[i]; 
  }

  return;     
}
//****************************************************************************80



void Xover ( int one, int two, int &seed ) {
  int i;
  int point;
  double t;
// 
//  Select the crossover point.
//
  point = i4_uniform_ab ( 0, NVARS - 1, seed );
//
//  Swap genes in positions 0 through POINT-1.
//
  for ( i = 0; i < point; i++ )
  {
    t                       = population[one].gene[i];
    population[one].gene[i] = population[two].gene[i];
    population[two].gene[i] = t;
  }

  return;
}
