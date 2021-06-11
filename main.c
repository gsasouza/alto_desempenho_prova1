#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>


int calculate_index(int i, int j, int k, int size_i, int size_j, int size_k) {
  return (size_j * size_k * i) + (size_k * j) + k;
}

int calculate_region_size(int size_city, int size_students) {
  return size_city * size_students;
}

int calculate_country_size(int size_region, int size_city, int size_students) {
  return size_region * size_city * size_students;
}

int compare_short(const void *a, const void *b) {
  return *((short int *) a) - *((short int *) b);
}

void sort_short(short int *start_pointer, int start, int size) {
  qsort(start_pointer + start, size, sizeof(short int), compare_short);
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

float calculate_mean_short(const short int *array, int start, int size) {
  int sum = 0;
  for (int i = start; i < start + size; i++) {
    sum += array[i];
  }
  return (float) sum / (float) size;
}

// standard deviation
float calculate_standard_deviation_short(short int *array, int start, int size, float mean) {
  float dp = (float) 0.0;
  for (int i = start; i < start + size; ++i) {
    dp += pow(array[i] - mean, 2);
  }
  return sqrt(dp / (float) size);
}

void calculate_city_stats(
  int size_regions,
  int size_cities,
  int size_students,
  short int *students_grade,
  float *city_median,
  float *city_mean,
  short int *city_highest,
  short int *city_lowest,
  float *city_standard_deviation
) {
  for (int i = 0; i < size_regions; i++) {
    for (int j = 0; j < size_cities; j++) {
      int position = calculate_index(0, i, j, 0, size_regions, size_cities);
      int city_start = calculate_index(i, j, 0, size_regions, size_cities, size_students);
      sort_short(students_grade, city_start, size_students); // sort each array slice (city)
      city_median[position] = calculate_median_short(students_grade, city_start, size_students);
      city_mean[position] = calculate_mean_short(students_grade, city_start, size_students);
      city_highest[position] = calculate_highest_short(students_grade, city_start, size_students);
      city_lowest[position] = calculate_lowest_short(students_grade, city_start);
      city_standard_deviation[position] = calculate_standard_deviation_short(students_grade, city_start,
                                                                      size_students, city_mean[j]);

    }
  }
}

void calculate_region_stats(
  int size_regions,
  int size_cities,
  int size_students,
  short int *students_grade,
  float *region_median,
  float *region_mean,
  short int *region_highest,
  short int *region_lowest,
  float *region_standard_deviation
) {

  for (int i = 0; i < size_regions; i++) {
    int region_start = calculate_index(i, 0, 0, size_regions, size_cities, size_students);
    int region_size = calculate_region_size(size_cities, size_students);
    sort_short(students_grade, region_start, region_size); // sort array slice (region)
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
  float *country_mean,
  short int *country_highest,
  short int *country_lowest,
  float *country_standard_deviation
) {
  int country_size = calculate_country_size(size_regions, size_cities, size_students);
  sort_short(students_grade, 0, country_size); // sort array slice (country)
  *country_median = calculate_median_short(students_grade, 0, country_size);
  *country_mean = calculate_mean_short(students_grade, 0, country_size);
  *country_highest = calculate_highest_short(students_grade, 0, country_size);
  *country_lowest = calculate_lowest_short(students_grade, 0);
  *country_standard_deviation = calculate_standard_deviation_short(students_grade, 0,
                                                                   country_size, *country_mean);
}

void pretty_print_output(
  int size_regions,
  int size_cities,
  int size_students,
  float *region_median,
  const float *region_mean,
  short int *region_highest,
  short int *region_lowest,
  float *region_standard_deviation,
  float *city_median,
  float *city_mean,
  short int *city_highest,
  short int *city_lowest,
  float *city_standard_deviation,
  float country_median,
  float country_mean,
  short int country_highest,
  short int country_lowest,
  float country_standard_deviation,
  double response_time
) {
  float best_region[2] = { -1, -1 }, best_city[3] = { -1, -1, -1}; // value = 0, region = 1, city = 2;
  for (int i = 0; i < size_regions; i++) {
    if (region_mean[i] > best_region[0]) {
      best_region[0] = region_mean[i];
      best_region[1] = (float) i;
    }
    for (int j = 0; j < size_cities; j++) {
      int position = calculate_index(0, i, j, 0,size_regions, size_cities);
      if (city_mean[position] > best_city[0]) {
        best_city[0] = city_mean[position];
        best_city[1] = (float) i;
        best_city[2] = (float) j;
      }

      printf("Reg %d - Cid %d: menor: %d, maior: %d, mediana: %.2f, média: %.2f e DP: %.2f\n", i, j,
             city_lowest[position],
             city_highest[position], city_median[position], city_mean[position], city_standard_deviation[position]);
    }
    printf("\n\n");
  }

  for(int i = 0; i < size_regions; i++) {
    printf("Reg %d: menor: %d, maior: %d, mediana: %.2f, média: %.2f e DP: %.2f\n", i, region_lowest[i], region_highest[i], region_median[i], region_mean[i], region_standard_deviation[i]);
  }
  printf("\n\n");

  printf("Brasil: menor: %d, maior: %d, mediana %.2f, média: %.2f e DP: %.2f \n\n", country_lowest, country_highest,
         country_median, country_mean, country_standard_deviation);

  printf("Melhor região: Região %0.f\n", best_region[1]);
  printf("Melhor cidade: Região %0.f, Cidade %0.f\n\n", best_city[1], best_city[2]);
  printf("Tempo de resposta sem considerar E/S, em segundos: %lf", response_time);

}


int main() {
  double response_time;
  int seed, size_students, size_regions, size_cities;
  scanf("%d %d %d %d", &size_regions, &size_cities, &size_students, &seed);
  short int *students_grade = malloc(sizeof(short int) * size_regions * size_cities * size_students);
  float *city_median = malloc(sizeof(float *) * size_regions * size_cities);
  float *city_mean = malloc(sizeof(float *) * size_regions * size_cities);
  short int *city_highest = malloc(sizeof(short int *) * size_regions * size_cities);
  short int *city_lowest = malloc(sizeof(short int *) * size_regions * size_cities);
  float *city_standard_deviation = malloc(sizeof(float *) * size_regions * size_cities);
  float *region_median = (float *) malloc(sizeof(float) * size_regions);
  float *region_mean = (float *) malloc(sizeof(float) * size_regions);
  short int *region_highest = (short int *) malloc(sizeof(short int) * size_regions);
  short int *region_lowest = (short int *) malloc(sizeof(short int) * size_regions);
  float *region_standard_deviation = (float *) malloc(sizeof(float) * size_regions);
  float country_median, country_mean, country_standard_deviation;
  short int country_highest, country_lowest;
  srand(seed);

  response_time = omp_get_wtime();

  for (int i = 0; i < size_regions; i++) {
    for (int j = 0; j < size_cities; j++) {
      for (int k = 0; k < size_students; k++) {
        students_grade[calculate_index(i, j, k, size_regions, size_cities, size_students)] = rand() % 100;
      }
    }
  }

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

  return 0;
}