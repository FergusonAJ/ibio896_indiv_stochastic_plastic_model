#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <random>
#include <fstream>

struct Predator{
  double x = 0;
  double y = 0;
  double a, b;
  std::array<double, 2> attack_rates;
  std::array<double, 2> handling_times;
  std::array<double, 2> functional_responses;
  std::array<double, 2> foraging_rates;
  std::array<double, 2> conv_efficiencies;
  std::array<double, 2> profitabilities;
  double prob_hab_1;
  double birth_val;
  Predator() { }
  Predator(double x_, double y_) : x(x_), y(y_){
    a = x - y;
    if(a < -1) a = -1;
    else if(a > 1) a = 1;
    b = x + y;
    if(b < -1) b = -1;
    else if(b > 1) b = 1;
  }
};

class Simulation{
  protected:
  // General variables
  size_t n_1;
  size_t n_2;
  std::vector<Predator> pred_vec;
  std::default_random_engine rand;
  // Constants
  const double growth_prey = 0.9;         // r_i = Intrinsic growth rate of both prey pops
  const double birth_prey = 1.0;          // b_i = Per-capita birth rate of both prey pops
  const double death_prey = 0.1;          // d_i = Per-capita death rate of both prey pops
  const double carrying_cap_prey = 2333;  // k_i = Carrying capacity of both prey pops
  const double max_attack_rate = 0.004;   // a_hat_i = Maximal attack rate of predator on prey i
  const double min_handling_time = 0.4;   // h_hat_i = Minimal handling time of prey i by predator
  const double max_conv_eff_prey_1 = 0.6; // c_hat_1 = Predator's maximal conversion efficeny of prey 1
  const double max_conv_eff_prey_2 = 0.2; // c_hat_2 = Predator's maximal conversion efficeny of prey 2
  const double x_mut_prob = 0.01;         // mu_x = Mutation probability for x (not mean)
  const double y_mut_prob = 0.01;         // mu_x = Mutation probability for y (not mean)
  const double x_mut_std = 0.1;           // sigma_x = Standard deviation of x mutation distribution
  const double y_mut_std = 0.1;           // sigma_y = Standard deviation of y mutation distribution
  const double a_std = 1;                 // sigma_a = Standard deviation of "a" distribution
  const double a_var = a_std * a_std;     // sigma^2_a = Variance of "a" distribution
  const double h_std = 1;                 // sigma_h = Standard deviation of "h" distribution
  const double h_var = h_std * h_std;     // sigma^2_h = Variance of a "h" disttribution
  
  public: 
  double death_pred = 0.22;         // m = Per-capita death rate for predators
  double plastic_cost = 0.1;        // k = Scaling factor for the cost of plasticity
  double habitat_sensitivity = 10;  // s = Habitat sensitivity
  size_t init_n_1 = 500;
  size_t init_n_2 = 500;
  size_t init_n_pred = 500;
  size_t cur_n_1;
  size_t cur_n_2;
  double time_elapsed = 0;
  double time_max = 10000;
  std::string file_prefix;

  protected:

  void UpdatePredator(Predator & pred){
    // Attack rates
    pred.attack_rates[0] = 
        max_attack_rate * exp(-1.0 * ( ((-1 - pred.a) * (-1 - pred.a)) / (2 * a_std) )); 
    pred.attack_rates[1] = 
        max_attack_rate * exp(-1.0 * ( ((1 - pred.b) * (1 - pred.b)) / (2 * a_std) )); 
    // Handling times
    pred.handling_times[0] = 
      min_handling_time + (1 - exp(-1.0 * ( ((-1 - pred.a) * (-1 - pred.a)) / (2 * h_std) ))); 
    pred.handling_times[1] = 
      min_handling_time + (1 - exp(-1.0 * ( ((1 - pred.b) * (1 - pred.b)) / (2 * h_std) ))); 
    // Functional responses
    pred.functional_responses[0] = 
      (pred.attack_rates[0] * cur_n_1) / (1 + pred.attack_rates[0] * pred.handling_times[0] * cur_n_1);  
    pred.functional_responses[1] = 
      (pred.attack_rates[1] * cur_n_2) / (1 + pred.attack_rates[1] * pred.handling_times[1] * cur_n_2);  
    // Conversion efficiencies
    pred.conv_efficiencies[0] = max_conv_eff_prey_1 * (1 - pred.y * plastic_cost);
    pred.conv_efficiencies[1] = max_conv_eff_prey_2 * (1 - pred.y * plastic_cost);
    // Profitabilities
    pred.profitabilities[0] = (pred.conv_efficiencies[0] * pred.functional_responses[0]) / death_pred;  
    pred.profitabilities[1] = (pred.conv_efficiencies[1] * pred.functional_responses[1]) / death_pred;  
    // Habitat selection probability
    pred.prob_hab_1 = 
        1 / (1 + exp(-1 * habitat_sensitivity * (pred.profitabilities[0] - pred.profitabilities[1])));
    // Foraging rate
    pred.foraging_rates[0] = pred.prob_hab_1 * pred.functional_responses[0];
    pred.foraging_rates[1] = (1 - pred.prob_hab_1) * pred.functional_responses[1];
  }

  void FillPredatorPop(size_t num_preds){
    pred_vec.resize(num_preds);
    double x, y;
    std::normal_distribution<double> x_dist(0, 0.05);
    std::normal_distribution<double> y_dist(0, 0.05);
    for(size_t i = 0; i < num_preds; ++i){
      x = x_dist(rand);
      y = y_dist(rand);
      while(y < 0)
        y = y_dist(rand);
      pred_vec[i] = Predator(x,y);
    }
  }

  void Initialize(){
    cur_n_1 = init_n_1;
    cur_n_2 = init_n_2;
    FillPredatorPop(init_n_pred); 
    time_elapsed = 0;
  }

  public:
  Simulation(size_t seed) : rand(seed){
  }

  void Run(){
    Initialize();
    std::ofstream fp_pop_size;
    fp_pop_size.open(file_prefix + "pop_sizes.csv");
    // Create variables that will be reused throughout main loop
    double total_birth_prey_1;
    double total_birth_prey_2;
    double total_birth_pred;
    double total_death_prey_1;
    double total_death_prey_2;
    double total_death_pred;
    double total_sum;
    double p;
    fp_pop_size << "name,count,time" << std::endl;
    size_t update = 0;
    while(time_elapsed < time_max){
      if(update % 100 == 0){
        fp_pop_size << "predator," << pred_vec.size() << "," << time_elapsed << std::endl
                    << "prey1," << cur_n_1 << "," << time_elapsed << std::endl 
                    << "prey2," << cur_n_2 << "," << time_elapsed << std::endl; 
      }
      update++;
      //std::cout << time_elapsed << std::endl;
      total_birth_pred = 0;
      total_death_prey_1 = (growth_prey * cur_n_1 * cur_n_1) / carrying_cap_prey;
      total_death_prey_2 = (growth_prey * cur_n_2 * cur_n_2) / carrying_cap_prey;
      for(size_t pred_idx = 0; pred_idx < pred_vec.size(); ++pred_idx){
        Predator & pred = pred_vec[pred_idx];
        UpdatePredator(pred);
        pred.birth_val = 
            pred.conv_efficiencies[0] * pred.foraging_rates[0] + 
            pred.conv_efficiencies[1] * pred.foraging_rates[1];
        total_birth_pred += pred.birth_val;
        total_death_prey_1 += pred.foraging_rates[0];
        total_death_prey_2 += pred.foraging_rates[1];
        //std::cout << pred.prob_hab_1 << "  " << pred.foraging_rates[0] << "  " << pred.foraging_rates[1] << std::endl;
        //std::cout << pred.profitabilities[0] << "  " << pred.profitabilities[1] << std::endl;
      }
      total_birth_prey_1 = growth_prey * cur_n_1; 
      total_birth_prey_2 = growth_prey * cur_n_2; 
      total_death_pred = death_pred * pred_vec.size(); 
      total_sum = total_birth_prey_1 + total_birth_prey_2 + total_birth_pred + 
          total_death_prey_1 + total_death_prey_2 + total_death_pred;
      // Normalize so sum = 1
      total_birth_prey_1 /= total_sum;
      total_birth_prey_2 /= total_sum;
      total_birth_pred /= total_sum;
      total_death_prey_1 /= total_sum;
      total_death_prey_2 /= total_sum;
      total_death_pred /= total_sum;
      // Let some time pass with lambda = total (total = E in paper)
      std::exponential_distribution<double> exp_dist(total_sum);
      time_elapsed += exp_dist(rand);
      // Generate random number between 0 and 1
      p = rand() / (double)rand.max();
      if(p < total_birth_prey_1)
        cur_n_1++;
      else if((p -= total_birth_prey_1)  < total_birth_prey_2)
        cur_n_2++;
      else if((p -= total_birth_prey_2) < total_death_prey_1)
        cur_n_1--;
      else if((p -= total_death_prey_1) < total_death_prey_2)
        cur_n_2--;
      else if((p -= total_death_prey_2) < total_birth_pred){
        double p2 = (rand() / (double)rand.max()) * total_birth_pred * total_sum;
        for(size_t pred_idx = 0; pred_idx < pred_vec.size(); ++pred_idx){
          if(p2 < pred_vec[pred_idx].birth_val){
            double x = pred_vec[pred_idx].x;
            double y = pred_vec[pred_idx].y;
            if(rand() / (double)rand.max() < x_mut_prob){
              std::normal_distribution<double> x_dist(pred_vec[pred_idx].x, x_mut_std);
              x = x_dist(rand);
            }
            if(rand() / (double)rand.max() < y_mut_prob){
              std::normal_distribution<double> y_dist(pred_vec[pred_idx].y, y_mut_std);
              y = y_dist(rand);
              while(y < 0)
                y = y_dist(rand);
            }
            pred_vec.emplace_back(x, y);
            break;
          }
          else{
            p2 -= pred_vec[pred_idx].birth_val;
          }
        }
      }
      else if((p -= total_birth_pred) < total_death_pred){
        size_t idx = rand() % pred_vec.size();
        pred_vec.erase(pred_vec.begin() + idx);
      }
    }
    fp_pop_size.close();
    std::ofstream fp_pred_genotypes;
    fp_pred_genotypes.open(file_prefix + "final_pred_genotypes.csv");
    fp_pred_genotypes << "idx, x, y" << std::endl; 
    for(size_t pred_idx = 0; pred_idx < pred_vec.size(); ++pred_idx){
      fp_pred_genotypes << pred_idx << "," << pred_vec[pred_idx].x << "," << pred_vec[pred_idx].y 
          << std::endl;
    }
    fp_pred_genotypes.close();
  } 

  void Reseed(size_t seed){
    rand.seed(seed);
  }

};

void PrintHelp(){
  std::cout << "Thanks for using the StochasticPlastic simulation!" << std::endl;
  std::cout << "Here are the possible command line arguments:" << std::endl;
  std::cout << "\t-h    Shows this help screen " << std::endl;
  std::cout << "\t-m <double>    Sets the predator mortality rate" << std::endl;
  std::cout << "\t-k <double>    Sets the cost of plasticity" << std::endl;
  std::cout << "\t-s <double>    Sets the predator habitat sensitivity" << std::endl;
  std::cout << "\t-p <int>    Sets the initial size of the predator population" << std::endl;
  std::cout << "\t-n1 <int>    Sets the initial size of the prey 1 population" << std::endl;
  std::cout << "\t-n2 <int>    Sets the initial size of the prey 2 population" << std::endl;
  std::cout << "\t-t <double>    Sets the time simulation will run" << std::endl;
  std::cout << "\t-r <int>    Sets the random seed" << std::endl;
  std::cout << "\t-f <string>    Sets the filename prefix. Should end in /" << std::endl;
}

int main(int argc, char * argv[]){
  // Setup simulation
  Simulation sim(time(NULL));
  // Interpret command line arguments
  bool skip = false;
  for(size_t arg_idx = 1; arg_idx < argc; arg_idx++){
    if(skip){
      skip = false;
      continue;
    }
    std::string arg_str(argv[arg_idx]);
    if(arg_str == "-h"){ 
      PrintHelp();
      exit(0);
    }
    else if(arg_str == "-m"){  // Predator mortality
      skip = true;
      if(arg_idx == argc - 1){
        std::cerr << "Error! Expected a double after -m!" << std::endl;
        exit(1);
      }
      double mortality = std::stod(argv[arg_idx+1]); 
      sim.death_pred  = mortality; 
    }
    else if(arg_str == "-k"){  // Plasticity cost
      skip = true;
      if(arg_idx == argc - 1){
        std::cerr << "Error! Expected a double after -k!" << std::endl;
        exit(1);
      }
      double cost = std::stod(argv[arg_idx+1]); 
      sim.plastic_cost  = cost; 
    }
    else if(arg_str == "-s"){  // Habitat sensitivity
      skip = true;
      if(arg_idx == argc - 1){
        std::cerr << "Error! Expected a double after -s!" << std::endl;
        exit(1);
      }
      double sensitivity = std::stod(argv[arg_idx+1]); 
      sim.habitat_sensitivity = sensitivity; 
    }
    else if(arg_str == "-p"){  // Initial predator count
      skip = true;
      if(arg_idx == argc - 1){
        std::cerr << "Error! Expected an int after -p!" << std::endl;
        exit(1);
      }
      int num_pred = std::stoi(argv[arg_idx+1]); 
      sim.init_n_pred = num_pred; 
    }
    else if(arg_str == "-n1"){ // Initial prey 1 count
      skip = true;
      if(arg_idx == argc - 1){
        std::cerr << "Error! Expected an int after -n1!" << std::endl;
        exit(1);
      }
      int num_1 = std::stoi(argv[arg_idx+1]); 
      sim.init_n_1 = num_1; 
    }
    else if(arg_str == "-n2"){ // Initial prey 2 count
      skip = true;
      if(arg_idx == argc - 1){
        std::cerr << "Error! Expected an int after -n2!" << std::endl;
        exit(1);
      }
      int num_2 = std::stoi(argv[arg_idx+1]); 
      sim.init_n_2 = num_2; 
    }
    else if(arg_str == "-t"){  // Maximum time 
      skip = true;
      if(arg_idx == argc - 1){
        std::cerr << "Error! Expected a double after -t!" << std::endl;
        exit(1);
      }
      double time = std::stod(argv[arg_idx+1]); 
      sim.time_max = time; 
    }
    else if(arg_str == "-r"){  // Random seed
      skip = true;
      if(arg_idx == argc - 1){
        std::cerr << "Error! Expected an int after -r!" << std::endl;
        exit(1);
      }
      int seed = std::stoi(argv[arg_idx+1]); 
      sim.Reseed(seed); 
    }
    else if(arg_str == "-f"){  // Random seed
      skip = true;
      if(arg_idx == argc - 1){
        std::cerr << "Error! Expected a string after -f!" << std::endl;
        exit(1);
      }
      std::string prefix(argv[arg_idx+1]); 
      sim.file_prefix = prefix; 
    }
    else{
      std::cout << "Unknown command line argument: " << arg_str << std::endl;
      exit(1);
    }
  }
  // Run
  sim.Run();
  return 0;
}
