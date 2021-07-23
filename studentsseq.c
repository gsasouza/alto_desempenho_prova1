/**********************************************************************/
/*Grupo 5 -> Questão 4 PB1 CAD 2021
/*Caio Augusto Duarte Basso 	- NUSP 10801173
/*Carla Nunes da Cruz       	- NUSP 8479343
/*Gabriel Garcia Lorencetti     	- NUSP 10691891
/*Gabriel Santos Souza - NUSP 11208176
/*Giovana Daniele da Silva    	- NUSP 10692224
/**********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <string.h>

#define RANGE 101
#define T 8


int calculate_index(int i, int j, int k, int size_i, int size_j, int size_k) {
  return (size_j * size_k * i) + (size_k * j) + k;
}

int calculate_region_size(int size_city, int size_students) {
  return size_city * size_students;
}

int calculate_country_size(int size_region, int size_city, int size_students) {
  return size_region * size_city * size_students;
}

void counting_sort(short int *array, int begin, int size) {
  short int *ordered = malloc(sizeof(short int) * size);
  int i;
  unsigned int counting[RANGE]; // [0, 100] is the range of the grades 

  memset(counting, 0, sizeof(counting));

  for (i = begin; i < begin + size; i++) {
    counting[array[i]] += 1;
  }

  for (int i = 1; i <= RANGE; i++) {
    counting[i] += counting[i - 1];
  }

  for (int i = begin; i < begin + size; i++) {
    ordered[counting[array[i]] - 1] = array[i];
    counting[array[i]]--;
  }

  for (int i = begin, j = 0; i < begin + size && j < size; i++, j++) {
    array[i] = ordered[j];
  }
  free(ordered);
}

short int calculate_lowest_short(short int *array, int start) {
  return array[start];
}

short int calculate_highest_short(short int *array, int start, int size) {
  return array[start + size - 1];
}

float calculate_median_short(short int *array, int start, int size) {
  int half = size / 2;
  if (size % 2) return array[start + half];
  return (float) (array[start + half - 1] + array[start + half]) / 2;
}

double calculate_mean_short(const short int *array, int start, int size) {
  double sum = 0.0;
  for (int i = start; i < start + size; i++) {
    sum = sum + array[i];
  }

  return sum / (double) size;
}

// sample standard deviation
double calculate_standard_deviation_short(short int *array, int start, int size, double mean) {
  double dp = 0.0;
  for (int i = start; i < start + size; i++) {
    dp += pow(array[i] - mean, 2);
  }
  return sqrt(dp / (double) (size - 1));
}

void calculate_city_stats(
  int size_regions,
  int size_cities,
  int size_students,
  short int *students_grade,
  float *city_median,
  double *city_mean,
  short int *city_highest,
  short int *city_lowest,
  double *city_standard_deviation
) {
  for (int i = 0; i < size_regions; i++) {
    #pragma omp parallel for num_threads(T)
    for (int j = 0; j < size_cities; j++) {
      int position = calculate_index(0, i, j, 0, size_regions, size_cities);
      int city_start = calculate_index(i, j, 0, size_regions, size_cities, size_students);
      counting_sort(students_grade, city_start, size_students); // sort each array slice (city)

      city_median[position] = calculate_median_short(students_grade, city_start, size_students);
      city_mean[position] = calculate_mean_short(students_grade, city_start, size_students);
      city_highest[position] = calculate_highest_short(students_grade, city_start, size_students);
      city_lowest[position] = calculate_lowest_short(students_grade, city_start);
      city_standard_deviation[position] = calculate_standard_deviation_short(students_grade, city_start,
                                                                             size_students, city_mean[position]);
    }
  }
}

void calculate_region_stats(
  int size_regions,
  int size_cities,
  int size_students,
  short int *students_grade,
  float *region_median,
  double *region_mean,
  short int *region_highest,
  short int *region_lowest,
  double *region_standard_deviation
) {

  #pragma omp parallel for num_threads(T) shared(students_grade, size_students, region_median, region_mean, region_highest, region_lowest, region_standard_deviation)
  for (int i = 0; i < size_regions; i++) {
    int region_start = calculate_index(i, 0, 0, size_regions, size_cities, size_students);
    int region_size = calculate_region_size(size_cities, size_students);

    counting_sort(students_grade, region_start, region_size); // sort array slice (region)

    region_median[i] = calculate_median_short(students_grade, region_start, region_size);
    region_mean[i] = calculate_mean_short(students_grade, region_start, region_size);
    region_highest[i] = calculate_highest_short(students_grade, region_start, region_size);
    region_lowest[i] = calculate_lowest_short(students_grade, region_start);
    region_standard_deviation[i] = calculate_standard_deviation_short(students_grade, region_start,
                                                                      region_size, region_mean[i]);
  }
}

void calculate_country_stats(
  int size_regions,
  int size_cities,
  int size_students,
  short int *students_grade,
  float *country_median,
  double *country_mean,
  short int *country_highest,
  short int *country_lowest,
  double *country_standard_deviation
) {

  int country_size = calculate_country_size(size_regions, size_cities, size_students);

  counting_sort(students_grade, 0, country_size);

  *country_mean = calculate_mean_short(students_grade, 0, country_size);

  *country_median = calculate_median_short(students_grade, 0, country_size);

  *country_highest = calculate_highest_short(students_grade, 0, country_size);
  *country_lowest = calculate_lowest_short(students_grade, 0);
  *country_standard_deviation = calculate_standard_deviation_short(students_grade, 0,
                                                                   country_size, *country_mean);
}

float rounded(float x) {
  float rounded = (int) (x * 100 + .5);

  return (float) rounded / 100;
}

void pretty_print_output(
  int size_regions,
  int size_cities,
  int size_students,
  float *region_median,
  const double *region_mean,
  short int *region_highest,
  short int *region_lowest,
  double *region_standard_deviation,
  float *city_median,
  double *city_mean,
  short int *city_highest,
  short int *city_lowest,
  double *city_standard_deviation,
  float country_median,
  double country_mean,
  short int country_highest,
  short int country_lowest,
  double country_standard_deviation,
  double response_time
) {

  float best_region[2] = {-1, -1}, best_city[3] = {-1, -1, -1}; // value = 0, region = 1, city = 2;
//  for (int i = 0; i < size_regions; i++) {
//    if (region_mean[i] > best_region[0]) {
//      best_region[0] = region_mean[i];
//      best_region[1] = (float) i;
//    }
//    for (int j = 0; j < size_cities; j++) {
//      int position = calculate_index(0, i, j, 0, size_regions, size_cities);
//      if (city_mean[position] > best_city[0]) {
//        best_city[0] = city_mean[position];
//        best_city[1] = (float) i;
//        best_city[2] = (float) j;
//      }
//
//      printf("Reg %d - Cid %d: menor: %d, maior: %d, mediana: %.2f, média: %.2lf e DP: %.2lf\n", i, j,
//             city_lowest[position],
//             city_highest[position], rounded(city_median[position]), rounded(city_mean[position]),
//             rounded(city_standard_deviation[position]));
//    }
//
//    printf("\n\n");
//  }

//  for (int i = 0; i < size_regions; i++) {
//    printf("Reg %d: menor: %d, maior: %d, mediana: %.2f, média: %.2lf e DP: %.2lf\n", i, region_lowest[i],
//           region_highest[i], rounded(region_median[i]), rounded(region_mean[i]),
//           rounded(region_standard_deviation[i]));
//  }
//
//  printf("\n\n");
//
//  printf("Brasil: menor: %d, maior: %d, mediana %.2f, média: %.2lf e DP: %.2lf \n\n", country_lowest, country_highest,
//         rounded(country_median), rounded(country_mean), rounded(country_standard_deviation));
//
//  printf("\nMelhor região: Região %0.f\n", best_region[1]);
//  printf("Melhor cidade: Região %0.f, Cidade %0.f\n\n", best_city[1], best_city[2]);
  printf("\nTempo de resposta sem considerar E/S, em segundos: %.3lfs\n", response_time);
}


int main(int argc, char *const argv[]) {
//  omp_set_nested(1);
  double response_time;
  int seed, size_students, size_regions, size_cities;

  scanf("%d %d %d %d", &size_regions, &size_cities, &size_students, &seed);

  short int *students_grade = malloc(sizeof(short int) * size_regions * size_cities * size_students);

  float *city_median = malloc(sizeof(float *) * size_regions * size_cities);
  double *city_mean = malloc(sizeof(double *) * size_regions * size_cities);
  double *city_standard_deviation = malloc(sizeof(double *) * size_regions * size_cities);
  short int *city_highest = malloc(sizeof(short int *) * size_regions * size_cities);
  short int *city_lowest = malloc(sizeof(short int *) * size_regions * size_cities);

  float *region_median = (float *) malloc(sizeof(float) * size_regions);
  double *region_mean = (double *) malloc(sizeof(double) * size_regions);
  double *region_standard_deviation = (double *) malloc(sizeof(double) * size_regions);
  short int *region_highest = (short int *) malloc(sizeof(short int) * size_regions);
  short int *region_lowest = (short int *) malloc(sizeof(short int) * size_regions);

  float country_median;
  double country_mean, country_standard_deviation;
  short int country_highest, country_lowest;

  srand(seed);



  for (int i = 0; i < size_regions; i++) {
    for (int j = 0; j < size_cities; j++) {
      for (int k = 0; k < size_students; k++) {
        students_grade[calculate_index(i, j, k, size_regions, size_cities, size_students)] = rand() % 101;
      }
    }
  }


  response_time = omp_get_wtime();
  /*short int students_grade[] = { 30, 40, 20, 80, 85, 10,
                     10, 20, 30, 40, 50, 60,
                     60, 50, 40, 30, 20, 10,
                     70, 55, 35, 80, 95, 27,
                     35, 45, 25, 85, 90, 15,
                     15, 25, 35, 45, 55, 65,
                     65, 55, 45, 35, 25, 15,
                     75, 60, 40, 85, 100, 32,
                     20, 30, 10, 70, 75, 0,
                     0,  10, 20, 30, 40, 50,
                     50, 40, 30, 20, 10, 0,
                     60, 45, 25, 70, 85, 17 };*/


  calculate_city_stats(
    size_regions,
    size_cities,
    size_students,
    students_grade,
    city_median,
    city_mean,
    city_highest,
    city_lowest,
    city_standard_deviation
  );


  calculate_region_stats(
    size_regions,
    size_cities,
    size_students,
    students_grade,
    region_median,
    region_mean,
    region_highest,
    region_lowest,
    region_standard_deviation
  );

  calculate_country_stats(
    size_regions,
    size_cities,
    size_students,
    students_grade,
    &country_median,
    &country_mean,
    &country_highest,
    &country_lowest,
    &country_standard_deviation
  );

  response_time = omp_get_wtime() - response_time;
//
  pretty_print_output(
    size_regions,
    size_cities,
    size_students,
    region_median,
    region_mean,
    region_highest,
    region_lowest,
    region_standard_deviation,
    city_median,
    city_mean,
    city_highest,
    city_lowest,
    city_standard_deviation,
    country_median,
    country_mean,
    country_highest,
    country_lowest,
    country_standard_deviation,
    response_time
  );

  free(students_grade);

  free(city_median);
  free(city_mean);
  free(city_standard_deviation);
  free(city_highest);
  free(city_lowest);

  free(region_median);
  free(region_mean);
  free(region_standard_deviation);
  free(region_highest);
  free(region_lowest);

  return 0;
}
