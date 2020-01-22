//
//  main.cpp
//  Assignment1
//
//  Created by Sultanmurat on 9/16/19.
//  Copyright © 2019 Sultanmurat. All rights reserved.
//

#include<iostream>
#include<fstream>
#include<iomanip>
#include<vector>
#include"/Users/zhazira/Desktop/Evolution/Assignment1/Assignment1/gnuplot_i.hpp"
#include<algorithm>
#include<ctime>
#include<random>       // std::default_random_engine
#include<chrono>       // std::chrono::system_clock
#include<stdlib.h>
#include<map>
#include<array>
#include<cmath>

using namespace std;

unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
std::random_device rd; // obtain a random number from hardware
std::mt19937 eng(rd()); // seed the generator
std::uniform_real_distribution<> distr01(0.0, 1.0);



double x, y;
struct Point{
    double x_coord;
    double y_coord;
};

// parameters:

const int n_population = 50;
double mutationProbability = 0.9;
double mutationProbability2 = 0.9;
double crossoverProbability = 1;
int number_attempts = 300000;  // number of maximum generations

// variables
vector<vector<int>> population;
vector<double> distances;
vector<Point> points;
size_t n_points;
double best_short_distances[10000000] = {0};
double best_long_distances[10000000] = {0};
double errors[10000000] = {0};
double distance(Point a, Point b);

double total_distance(vector<int> cities);

bool compareFitness(array<double, 2> fitness1, array<double, 2> fitness2);
bool compareIndividuals(vector<int> ind1, vector<int> ind2);

vector<int> order_crossover(vector<int> parent1, vector<int> parent2, int cutPoint1, int cutPoint2);
vector<int> pmx_crossover(vector<int> parent1, vector<int> parent2, int cutPoint1, int cutPoint2);
vector<int> crossover_ends(vector<int> parent1, vector<int> parent2, int cutPoint1, int cutPoint2);

void wait_for_key();
void print_distances();
void circle_create(double radius, int n);

vector<int> mutateBig4(vector<int> child1V);
vector<int> mutateBig1(vector<int> child1V);
vector<int> mutateBig2(vector<int> child1V);
vector<int> mutateBig3(vector<int> child1V);
vector<int> mutateSmall1(vector<int> child1V);
vector<int> mutateSmall2(vector<int> child1V);
vector<int> mutateSmall3(vector<int> child1V);
vector<int> mutateEveryElement(vector<int> child1V);

int main(){
    circle_create(1, 12);
  
    // Reading points
    ifstream input;
    input.open("/Users/zhazira/Desktop/Evolution/Assignment1/Assignment1/circle.txt");
    
    while (input >> x >> y){
        Point a;
        a.x_coord = x;
        a.y_coord = y;
        points.push_back(a);
    }
    
    n_points = points.size();
    
    vector<int> cities;
    for (int i = 0; i < n_points; i++){
        cities.push_back(i);
    }
 
    // Outputting values
    ofstream circle_answer;
    circle_answer.open("/Users/zhazira/Desktop/Evolution/Assignment1/Assignment1/circle_answer.txt");
    
    
    
    
    // Creating initial random population
    
    for (int i = 0; i < n_population; i++){
        shuffle(cities.begin(), cities.end(), default_random_engine(seed)); //  -- Random
        population.push_back(cities);
    }
    int j = 0;
    
    // random search
//    vector<vector<Point>> evaluations;
//    int max = 0;
//    int min = 0;
//    int index_min = 0;
//    int index_max = 0;
//    while (j <= number_attempts){
//        shuffle (points.begin(), points.end(), default_random_engine(seed));
//        vector<Point> temp_points = points;
//        evaluations.push_back(temp_points);
//        double total_distance = 0;
//        for (int i = 0; i < n_points-1; i++) {
//            double current_distance = distance(points[i], points[i+1]);
//            total_distance += current_distance;
//        }
//        distances.push_back(total_distance);
//        //if (total_distance < 500) break;
//        if (j==0 || total_distance < min){
//            min = total_distance;
//            index_min = j;
//        }
//        if (j==0 || total_distance > max){
//            max = total_distance;
//            index_max = j;
//        }
//        best_short_distances[j] = min;
//        best_long_distances[j] = max;
//        j++;
//        cout << j << endl;
//    }
//    ofstream citiesShort;
//    ofstream citiesLong;
//    ofstream learningCurveShort;
//    ofstream learningCurveLong;
//    citiesShort.open("/Users/zhazira/Desktop/Evolution/Assignment1/Assignment1/citiesShort.txt");
//    citiesLong.open("/Users/zhazira/Desktop/Evolution/Assignment1/Assignment1/citiesLong.txt");
//    learningCurveShort.open("/Users/zhazira/Desktop/Evolution/Assignment1/Assignment1/curve_short.txt");
//    learningCurveLong.open("/Users/zhazira/Desktop/Evolution/Assignment1/Assignment1/curve_long.txt");
//
//    for (int i = 0; i < n_points; i++){
//        citiesShort << setprecision(12) << evaluations[index_min][i].x_coord << " " << setprecision(12) << evaluations[index_min][i].y_coord << endl;
//         citiesLong << setprecision(12) << evaluations[index_max][i].x_coord << " " << setprecision(12) << evaluations[index_max][i].y_coord << endl;
//    }
//    for (int i = 0; i < number_attempts; i++){
//        //output2 << i+1 << " " << distances[i] << endl;
//        learningCurveShort << i+1 << " " << best_short_distances[i] << endl;
//        learningCurveLong << i+1 << " " << best_long_distances[i] << endl;
//    }
//
//    input.close();
//    citiesShort.close();
//    citiesLong.close();
//    learningCurveLong.close();
//    learningCurveShort.close();
//
//    Gnuplot gp("lines");
//    gp.plotfile_xy("/Users/zhazira/Desktop/Evolution/Assignment1/Assignment1/citiesShort.txt", 1, 2, "Cities Short");
//    gp.showonscreen();
//
//    Gnuplot gp2("lines");
//    gp2.plotfile_xy("/Users/zhazira/Desktop/Evolution/Assignment1/Assignment1/citiesLong.txt", 1, 2, "Cities Long");
//    gp2.showonscreen();
//
//
//
//    Gnuplot gp3("lines");
//    gp3.set_ylabel("Fitness");
//    gp3.set_xlabel("Evalatuations");
//    gp3.set_title("Shortest and longest path learning curves for random search");
//    gp3.plotfile_xy("/Users/zhazira/Desktop/Evolution/Assignment1/Assignment1/curve_short.txt", 1, 2, "Learning curve shortest path");
//    gp3.plotfile_xy("/Users/zhazira/Desktop/Evolution/Assignment1/Assignment1/curve_long.txt", 1, 2, "Learning curve longest path");
//    gp3.showonscreen();


    
    //hill climber for shortest
//    while (j < number_attempts ){
//        double current_distances[n_population];
//        double total_distances = 0;
//        for (int i = 0; i < n_population; i++){
//            int position1 = (int) (distr01(eng)*n_points);
//            int position2 = (int) (distr01(eng)*n_points);
//            vector<int> temp_cities;
//            temp_cities = cities;
//            swap(temp_cities[position1], temp_cities[position2]);
//            double prev_distance = total_distance(cities);
//            current_distances[i] = prev_distance;
//            double new_distance =total_distance(temp_cities);
//            if (new_distance < prev_distance){
//                cities = temp_cities;
//                current_distances[i] = new_distance;
//            }
//
//            total_distances += current_distances[i];
//        }
//        double mean = total_distances/n_population;
//        double error_deviation = 0;
//        for (int i = 0; i < n_population; i++){
//            error_deviation = error_deviation + sqrt((current_distances[i]-mean)*(current_distances[i]-mean));
//        }
//        error_deviation = error_deviation/n_population;
//        errors[j] = error_deviation;
//        best_short_distances[j] = mean;
//        cout << j << "   Shortest:  " << setprecision(4) << best_short_distances[j] << endl;
//        j++;
//    }
//    cout << "The shortest distance using hill climber is: " << total_distance(cities) << endl;
    
//    //hill climber for longest
//     while (j < number_attempts ){
//         double current_distances[n_population];
//         double total_distances = 0;
//         for (int i = 0; i < n_population; i++){
//             int position1 = (int) (distr01(eng)*n_points);
//             int position2 = (int) (distr01(eng)*n_points);
//             vector<int> temp_cities;
//             temp_cities = cities;
//             swap(temp_cities[position1], temp_cities[position2]);
//             double prev_distance = total_distance(cities);
//             current_distances[i] = prev_distance;
//             double new_distance =total_distance(temp_cities);
//             if (new_distance > prev_distance){
//                 cities = temp_cities;
//                 current_distances[i] = new_distance;
//             }
//
//             total_distances += current_distances[i];
//         }
//         double mean = total_distances/n_population;
//         double error_deviation = 0;
//         for (int i = 0; i < n_population; i++){
//             error_deviation = error_deviation + sqrt((current_distances[i]-mean)*(current_distances[i]-mean));
//         }
//         error_deviation = error_deviation/n_population;
//         errors[j] = error_deviation;
//         best_short_distances[j] = mean;
//         cout << j << "   The longest:  " << setprecision(4) << best_short_distances[j] << endl;
//         j++;
//     }
//     cout << "The longest distance using hill climber is: " << total_distance(cities) << endl;


        
        
    int cutPoint1 = (int)n_points/4; //n_points/4 + distr01(eng)*n_points/4;
    int cutPoint2 = (int)(3*n_points/4); //n_points/2 + distr01(eng)*n_points/4;
     //Genetic algorithm, j - generation number
    double min;
    while (j < number_attempts){
        vector<vector<int>> next_population = population;

        // Get fitness values and sort them
        array<double, 2> fitness[n_population];
        double total_fitness = 0;
        for (int i = 0; i < n_population; i++){
            double current_total_distance = total_distance(population[i]);
            fitness[i][0] = current_total_distance;
            fitness[i][1] = i;
            total_fitness += fitness[i][0];
        }
        double mean = total_fitness/n_population;
        double error_deviation = 0;

        for (int i = 0; i < n_population; i++){
            error_deviation = error_deviation + sqrt((fitness[i][0]-mean)*(fitness[i][0]-mean));
        }

        error_deviation = error_deviation/n_population;
        //cout << "Error: " << error_deviation << endl;
        errors[j] = error_deviation;
        //sort(fitness, fitness + n_population, compareFitness);

        // Select parents for crossover using roulette selection


        double accumulated_fitness[n_population] = {0};
        accumulated_fitness[0] = fitness[0][0];
        for (int i = 1; i < n_population; i++){
            accumulated_fitness[i] = accumulated_fitness[i-1] + fitness[i][0];
        }
        //for (int i =0; i < n_population; i++){
        //    cout << fitness[i][0] << " " << fitness[i][1] << endl;
        //}



        std::uniform_real_distribution<> distr(0, total_fitness); // define the range
        std::uniform_int_distribution<> distr5(0, n_population-1);

        double p = total_fitness/n_population;
        int start = (int)(distr01(eng)*p);

        // Producing new children
        for (int m =0; m < n_population; m++){
            // Select using tournament selection
//            int k = n_population/4;
//            int r1 = 0;
//            int r2 = 0;
//            for (int i =0; i < k; i++){
//                int ind = (int) (distr01(eng)*n_population);
//                if (r1 == 0 || fitness[ind][0] < fitness[r1][0]){
//                    r1 = ind;
//                }
//            }
//            for (int i =0; i < k; i++){
//                int ind = (int) (distr01(eng)*n_population);
//                if (r2 == 0 || fitness[ind][0] < fitness[r2][0]){
//                    r2 = ind;
//                }
//            }
//            r1 = fitness[r1][1];
//            r2 = fitness[r2][1];

           // Stochastic selection
            int r1 = start+m*p/2;
            int r2 = start+m*p;
            if (p==0) {
                r1 = start/2;
                r2 = start;}
            // Selection roulette
            //int r1 = distr(eng);
            //int r2 = distr(eng);

            for (int i = 0; i < n_population; i++){
                if (accumulated_fitness[i] >= r1){
                    r1 = fitness[i][1];
                    break;
                }
            }
            for (int i = 0; i < n_population; i++){
                if (accumulated_fitness[i] >= r2){
                    r2 = fitness[i][1];
                    break;
                }
            }

            // making crossover between parent1 and parent2 with PMX
            vector<int> parent1 = population[r1];
            vector<int> parent2 = population[r2];
            //const size_t n = n_points;




            // pushing children values into new population
            vector<int> child1V;
//            vector<int> child2V;
            child1V = order_crossover(parent1, parent2, cutPoint1, cutPoint2);
            //child1V = crossover_ends(parent1, parent2, cutPoint1, cutPoint2);
            // child2V = crossover_ends(parent2, parent1, cutPoint1, cutPoint2);
             //child1V = pmx_crossover(parent1, parent2, cutPoint1, cutPoint2);
           // child1V = mutateBig1(child1V);
            //child1V = mutateBig2(child1V);
            child1V = mutateBig3(child1V);
            //child1V = mutateBig4(child1V);
            //child1V = mutateSmall1(child1V);
            //child1V = mutateSmall2(child1V);
            //child1V = mutateSmall3(child1V);
//            child2V = mutateBig1(child2V);
//            child2V = mutateBig2(child2V);
//            child2V = mutateBig3(child2V);
//            child2V = mutateBig4(child2V);
//            child2V = mutateSmall1(child2V);
//            child2V = mutateSmall2(child2V);
//            child1V = mutateSmall3(child1V);
//            child2V = mutateSmall3(child2V);

            next_population.push_back(child1V);
            //next_population.push_back(child2V);

        }

        sort(next_population.begin(), next_population.end(), compareIndividuals);

        for (int i = 0; i < n_population; i++){
            next_population.pop_back();
        }

        population = next_population;

        if (total_distance(population[0]) < min || j==0){
            min = total_distance(population[0]);
        }
        best_short_distances[j] = min;
        double diff = total_distance(population[n_population-1])-total_distance(population[0]);
        cout << "Gen: " << j << "   Min: " << setprecision(3) << min << "   Diff: " << setprecision(5) << diff
         <<endl;
        j++;
        if (diff == 0) break;


    }
    print_distances();
    //cout << " Number of points" << n_points;
    Gnuplot gp("linespoints");
    for (int i = 0; i < n_points; i++){
          //output for genetic algorithm
          circle_answer << points[population[0][i]].x_coord << " " << points[population[0][i]].y_coord << endl ;
            //output for hill climber
            //circle_answer << points[cities[i]].x_coord << " " << points[cities[i]].y_coord << endl;
    }
    circle_answer << points[population[0][0]].x_coord << " " << points[population[0][0]].y_coord << endl;
    //circle_answer << points[cities[0]].x_coord << " " << points[cities[0]].y_coord << endl;
    gp.set_title("Shortest path found using GA with ordered crossover");
    gp.plotfile_xy("/Users/zhazira/Desktop/Evolution/Assignment1/Assignment1/circle_answer.txt", 1, 2, "Path");
    gp.showonscreen();

    ofstream learningCurveShort;
    learningCurveShort.open("/Users/zhazira/Desktop/Evolution/Assignment1/Assignment1/curve_short.txt");
       for (int i = 0; i < j; i++){
              //output2 << i+1 << " " << distances[i] << endl;
              learningCurveShort << i+1 << " " << best_short_distances[i] << " "  << errors[i] << endl;
              //learningCurveLong << i+1 << " " << best_long_distances[i] << endl;
          }
        learningCurveShort.close();

    Gnuplot gp3("lines");
    gp3.set_xlabel("Evaluations");
    gp3.set_ylabel("Fitness");
    gp3.set_title("Learning curve of GA for the shortest path with ordered crossover");
    gp3.plotfile_xy_err("/Users/zhazira/Desktop/Evolution/Assignment1/Assignment1/curve_short.txt", 1, 2, 3, "Learning curve with GA for the shortest Path with ordered crossover");

    //gp3.plotfile_xy("/Users/zhazira/Desktop/Evolution/Assignment1/Assignment1/curve_long.txt", 1, 2, "Learning curve long");
    gp3.showonscreen();

    wait_for_key();
    return 0;
}
double distance(Point a, Point b){
    double res = sqrt((a.x_coord-b.x_coord)*(a.x_coord-b.x_coord)+(a.y_coord-b.y_coord)*(a.y_coord-b.y_coord));
    return res;
}

double total_distance(vector<int> cities){
    double res = 0;
    const size_t n_points = points.size();
    
    for (int i = 0; i < n_points-1; i++){
        res += distance(points[cities[i]], points[cities[i+1]]);
    }
    res+= distance(points[cities[0]], points[cities[n_points-1]]);
    return res;
}

bool compareFitness(array<double, 2> fitness1, array<double, 2> fitness2)
{
    return (fitness1[0] < fitness2[0]);
}

bool compareIndividuals(vector<int> ind1, vector<int> ind2){
    return total_distance(ind1) < total_distance(ind2);
}

vector<int> order_crossover(vector<int> parent1, vector<int> parent2, int cutPoint1, int cutPoint2){
    vector<int> child1V;
    if (distr01(eng) >= crossoverProbability){
        int n = parent1.size();
        int child1[n];
        bool inChild1[n];
        for (int i = 0; i < n; i++){
            inChild1[i] = false;
            child1[i] = 0;
        }
        for (int i = cutPoint1; i < cutPoint2; i++){
            child1[i] = parent1[cutPoint2 - i-1];
            inChild1[child1[i]] = true;
        }
        int j = cutPoint2;
        int k = cutPoint2;
        while (k != cutPoint1 ){
            if (! inChild1[parent2[j]]){
                child1[k] = parent2[j];
                k++;
                if (k == n){
                    k = 0;
                }
            }
            j++;
            if (j == n){
                j = 0;
            }
        }
        for (int i = 0; i < n; i++){
            child1V.push_back(child1[i]);
        }
        return child1V;
    } else {
        return parent1;
    }
    
}

vector<int> pmx_crossover(vector<int> parent1, vector<int> parent2, int cutPoint1, int cutPoint2){
    if (distr01(eng) <= crossoverProbability){
        int n = parent1.size();
        map<int, int> pmx;
        for (int i = cutPoint1; i < cutPoint2; i++){
            pmx[parent1[i]] = parent2[i];
        }
        // producing Child1
        int child1[n];
        bool inChild1[n];
        for (int i = 0; i < n; i++){
            inChild1[i] = false;
            child1[i] = 0;
        }
        
        for (int i = cutPoint1; i < cutPoint2; i++){
            child1[i] = parent1[i];
            inChild1[child1[i]] = true;
        }
        
        for (int i = cutPoint2; i < n; i++){
            bool isInChild = inChild1[parent2[i]];
            if (!isInChild){
                child1[i] = parent2[i];
                inChild1[child1[i]] = true;
            }
            else {
                int x = pmx[parent2[i]];
                while (inChild1[x]==true){
                    x = pmx[x];
                }
                child1[i] = x;
                inChild1[x] = true;
            }
        }
        for (int i = 0; i < cutPoint1; i++){
            bool isInChild = inChild1[parent2[i]];
            if (!isInChild){
                child1[i] = parent2[i];
                inChild1[child1[i]] = true;
            }
            else {
                int x = pmx[parent2[i]];
                while (inChild1[x]==true){
                    x = pmx[x];
                }
                child1[i] = x;
                inChild1[x] = true;
            }
        }
        vector<int> child1V;
        
        for (int i = 0; i < n; i++){
            child1V.push_back(child1[i]);
        }
        return child1V;
    }
    else {
        return parent1;
    }
}

vector<int> crossover_ends(vector<int> parent1, vector<int> parent2, int cutPoint1, int cutPoint2){
    if (distr01(eng) <= crossoverProbability){
        int n = parent1.size();
        map<int, int> pmx;
        for (int i = cutPoint1; i < cutPoint2; i++){
            pmx[parent1[i]] = parent2[i];
        }
        // producing Child1
        int child1[n];
        bool inChild1[n];
        for (int i = 0; i < n; i++){
            inChild1[i] = false;
            child1[i] = 0;
        }
        
        for (int i = 0; i < cutPoint1; i++){
            child1[i] = parent1[i];
            inChild1[child1[i]] = true;
        }
        for (int i = cutPoint2; i < n; i++){
            child1[i] = parent1[i];
            inChild1[child1[i]] = true;
        }
        
        int j = 0;
        int k = cutPoint1;
        while (k != cutPoint2 ){
            if (! inChild1[parent2[j]]){
                child1[k] = parent2[j];
                k++;
            }
            j++;
        }
        vector<int> child1V;
        for (int i = 0; i < n; i++){
            child1V.push_back(child1[i]);
        }
        return child1V;
    }
    else {
        return parent1;
    }
}

void print_distances(){
    for (int i = 0; i < n_population; i++){
        cout << total_distance(population[i]) << endl;
    }
}

void circle_create(double radius, int n){
    double x[n];
    double y[n];
    ofstream output_circle;
    output_circle.open("/Users/zhazira/Desktop/Evolution/Assignment1/Assignment1/circle.txt");
    for (int i = 0; i < n-1; i++){
        x[i] = distr01(eng)*radius;
        y[i] = sqrt(radius*radius-x[i]*x[i]);
        output_circle << x[i] << " " << y[i] << endl;
        output_circle << -x[i] << " " << -y[i] << endl;
        output_circle << -x[i] << " " << y[i] << endl;
        output_circle << x[i] << " " << -y[i] << endl;
    }
    output_circle << -radius << " " << 0 << endl;
    output_circle << radius << " " << 0 << endl;
    output_circle << 0 << " " << -radius << endl;
       output_circle << 0 << " " << radius << endl;
}

vector<int> mutateBig4(vector<int> child1V){
    int n_points = child1V.size();
    
    for (int i =0; i < n_points; i++){
        if (distr01(eng) <= mutationProbability2){
            int position1 = i;
            int position2 = (int) (distr01(eng)*n_points);
            std::swap(child1V[position1], child1V[position2]);
        }
    }
    return child1V;
}

vector<int> mutateBig1(vector<int> child1V){
    int n_points = child1V.size();
    vector<int> temp_child = child1V;
    if (distr01(eng) <= mutationProbability){
        int c = (int)(distr01(eng)*n_points/2);
        for (int i =0; i < c; i++){
            std::swap(child1V[i], child1V[n_points-c+i]);
        }
    }
    return child1V;
}

vector<int> mutateBig2(vector<int> child1V){
    int n_points = child1V.size();
    vector<int> temp_child = child1V;
    if (distr01(eng) < mutationProbability){
        int c = (int)(distr01(eng)*n_points/4);
        for (int i =0; i < c; i++){
            std::swap(child1V[i], child1V[c+i]);
        }
    }
    return child1V;
}

vector<int> mutateBig3(vector<int> child1V){
    int n_points = child1V.size();
    vector<int> temp_child = child1V;
    if (distr01(eng) < mutationProbability){
        int c1 = (int)(distr01(eng)*n_points/2);
        int c2 = (int)(distr01(eng)*n_points/2+n_points/2);
        for (int i =0; i < int((c2-c1)/2); i++){
            std::swap(child1V[c1+i], child1V[c2-i]);
        }
    }
    return child1V;
}

vector<int> mutateSmall1(vector<int> child1V){
    int n_points = child1V.size();
    vector<int> temp_child = child1V;
    if (distr01(eng) <= mutationProbability){
        int c = (int)(distr01(eng)*n_points/2);
        for (int i =0; i < c; i++){
            std::swap(child1V[i], child1V[n_points-i-1]);
        }
    }
    return child1V;
}

vector<int> mutateSmall2(vector<int> child1V){
    int n_points = child1V.size();
    vector<int> temp_child = child1V;
    if (distr01(eng) < mutationProbability){
        int c = (int)(distr01(eng)*n_points-1);
        for (int i =0; i < (int)c/2; i++){
            std::swap(child1V[i], child1V[c-1-i]);
        }
    }
    return child1V;
}

vector<int> mutateSmall3(vector<int> child1V){
    int n_points = child1V.size();
    vector<int> temp_child = child1V;
    if (distr01(eng) < mutationProbability){
        int c = (int)(distr01(eng)*n_points-1);
        std::swap(child1V[c], child1V[c+1]);
    }
    return child1V;
}


void wait_for_key ()
{
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__TOS_WIN__)  // every keypress registered, also arrow keys
    cout << endl << "Press any key to continue..." << endl;
    
    FlushConsoleInputBuffer(GetStdHandle(STD_INPUT_HANDLE));
    _getch();
#elif defined(unix) || defined(__unix) || defined(__unix__) || defined(__APPLE__)
    cout << endl << "Press ENTER to continue..." << endl;
    
    std::cin.clear();
    std::cin.ignore(std::cin.rdbuf()->in_avail());
    std::cin.get();
#endif
    return;
}

