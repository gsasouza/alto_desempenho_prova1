#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <string.h>
#include <mpi.h>

#define RANGE 101
#define T 16

int calculate_index(int i, int j, int k, int size_i, int size_j, int size_k) {
  return (size_j * size_k * i) + (size_k * j) + k;
}

int calculate_region_size(int size_city, int size_students) {
  return size_city * size_students;
}

int calculate_country_size(int size_region, int size_city, int size_students) {
  return size_region * size_city * size_students;
}

void counting_sort(int *array, int begin, int size) {
  int *ordered = malloc(sizeof(int) * size);
  unsigned int counting[RANGE]; // [0, 100] is the range of the grades 

  memset(counting, 0, sizeof(counting));

  for (int i = begin; i < begin + size; i++) {
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

int calculate_lowest_short(int *array, int start) {
  return array[start];
}

int calculate_highest_short(int *array, int start, int size) {
  return array[start + size - 1];
}

float calculate_median_short(int *array, int start, int size) {
  int half = size / 2;
  if (size % 2) return array[start + half];
  return (float) (array[start + half - 1] + array[start + half]) / 2;
}

double calculate_mean_short(const int *array, int start, int size) {
  double sum = 0.0;
  for (int i = start; i < start + size; i++) {
    sum = sum + array[i];
  }

  return sum / (double) size;
}

// sample standard deviation
double calculate_standard_deviation_short(int *array, int start, int size, double mean) {
  double dp = 0.0;

  for (int i = start; i < start + size; i++) {
    dp = dp + pow(array[i] - mean, 2);
  }
  return sqrt(dp / (double) (size - 1));
}

void calculate_city_stats(
  int size_regions,
  int size_cities,
  int size_students,
  int *students_grade,
  float *city_median,
  double *city_mean,
  int *city_highest,
  int *city_lowest,
  double *city_standard_deviation
) {
#pragma omp parallel for num_threads(T)
  for (int i = 0; i < size_regions; i++) {
    for (int j = 0; j < size_cities; j++) {
      int position = calculate_index(0, i, j, 0, size_regions, size_cities);
      int city_start = calculate_index(i, j, 0, size_regions, size_cities, size_students);
      counting_sort(students_grade, city_start, size_students); // sort each array slice (city)
      city_median[position] = calculate_median_short(students_grade, city_start, size_students);
      city_highest[position] = calculate_highest_short(students_grade, city_start, size_students);
      city_lowest[position] = calculate_lowest_short(students_grade, city_start);

      city_mean[position] = calculate_mean_short(students_grade, city_start, size_students);

      city_standard_deviation[position] = calculate_standard_deviation_short(students_grade, city_start,
                                                                             size_students, city_mean[position]);

    }
  }
}

void calculate_region_stats(
  int size_regions,
  int size_cities,
  int size_students,
  int *students_grade,
  float *region_median,
  double *region_mean,
  int *region_highest,
  int *region_lowest,
  double *region_standard_deviation
) {
#pragma omp parallel for num_threads(T)
  for (int i = 0; i < size_regions; i++) {
    int region_start = calculate_index(i, 0, 0, size_regions, size_cities, size_students);
    int region_size = calculate_region_size(size_cities, size_students);
    counting_sort(students_grade, region_start, region_size); // sort array slice (region)
    region_median[i] = calculate_median_short(students_grade, region_start, region_size);
    region_highest[i] = calculate_highest_short(students_grade, region_start, region_size);
    region_lowest[i] = calculate_lowest_short(students_grade, region_start);
    region_mean[i] = calculate_mean_short(students_grade, region_start, region_size);
    region_standard_deviation[i] = calculate_standard_deviation_short(students_grade, region_start,
                                                                      region_size, region_mean[i]);

  }
}

void calculate_country_stats(
  int size_regions,
  int size_cities,
  int size_students,
  int *students_grade,
  float *country_median,
  double *country_mean,
  int *country_highest,
  int *country_lowest,
  double *country_standard_deviation
) {

  int country_size = calculate_country_size(size_regions, size_cities, size_students);
#pragma omp parallel num_threads(T) firstprivate(students_grade) shared(country_size, country_median, country_mean, country_highest, country_lowest, country_standard_deviation)
  {
#pragma omp single
    {
#pragma omp task
      {
        counting_sort(students_grade, 0, country_size);
        *country_median = calculate_median_short(students_grade, 0, country_size);
        *country_highest = calculate_highest_short(students_grade, 0, country_size);
        *country_lowest = calculate_lowest_short(students_grade, 0);
      }
#pragma omp task
      {
        *country_mean = calculate_mean_short(students_grade, 0, country_size);
        *country_standard_deviation = calculate_standard_deviation_short(students_grade, 0,
                                                                         country_size, *country_mean);
      };
    }
  }

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
  int *region_highest,
  int *region_lowest,
  double *region_standard_deviation,
  float *city_median,
  double *city_mean,
  int *city_highest,
  int *city_lowest,
  double *city_standard_deviation,
  float country_median,
  double country_mean,
  int country_highest,
  int country_lowest,
  double country_standard_deviation,
  double response_time
) {

  float best_region[2] = {-1, -1}, best_city[3] = {-1, -1, -1}; // value = 0, region = 1, city = 2;
  for (int i = 0; i < size_regions; i++) {
    if (region_mean[i] > best_region[0]) {
      best_region[0] = region_mean[i];
      best_region[1] = (float) i;
    }
    for (int j = 0; j < size_cities; j++) {
      int position = calculate_index(0, i, j, 0, size_regions, size_cities);
      if (city_mean[position] > best_city[0]) {
        best_city[0] = city_mean[position];
        best_city[1] = (float) i;
        best_city[2] = (float) j;
      }

      printf("Reg %d - Cid %d: menor: %d, maior: %d, mediana: %.2f, média: %.2lf e DP: %.2lf\n", i, j,
             city_lowest[position],
             city_highest[position], rounded(city_median[position]), rounded(city_mean[position]),
             rounded(city_standard_deviation[position]));
    }

    printf("\n\n");
  }

  for (int i = 0; i < size_regions; i++) {
    printf("Reg %d: menor: %d, maior: %d, mediana: %.2f, média: %.2lf e DP: %.2lf\n", i, region_lowest[i],
           region_highest[i], rounded(region_median[i]), rounded(region_mean[i]),
           rounded(region_standard_deviation[i]));
  }

  printf("\n\n");

  printf("Brasil: menor: %d, maior: %d, mediana %.2f, média: %.2lf e DP: %.2lf \n\n", country_lowest, country_highest,
         rounded(country_median), rounded(country_mean), rounded(country_standard_deviation));

  printf("\nMelhor região: Região %0.f\n", best_region[1]);
  printf("Melhor cidade: Região %0.f, Cidade %0.f\n\n", best_city[1], best_city[2]);
  printf("\nTempo de resposta sem considerar E/S, em segundos: %.3lfs\n", response_time);
}


void calculate_counts_displs(int *send_counts, int *displs, int *count_subset, int comm_size, int total_count,
                             int row_size) {
  int worker_processes = (comm_size - 1);
  int remainder = total_count % worker_processes;
  int count_per_process = total_count / worker_processes;
  int size_per_process = row_size * count_per_process;
  int sum = 0;
  displs[0] = 0;
  send_counts[0] = 0;
  count_subset[0] = 0;
  for (int i = 1; i < comm_size; i++) {
    send_counts[i] = size_per_process;
    count_subset[i] = count_per_process;
    if (remainder > 0) {
      send_counts[i] += row_size;
      count_subset[i] += 1;
      remainder--;
    }
    displs[i] = sum;
    sum += send_counts[i];
  }
}

int main(int argc, char *const argv[]) {
  MPI_Init(NULL, NULL);

  int rank, comm_size, chunk_size;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

  double response_time;
  int seed, size_students, size_regions, size_cities;

  // variables containing entire data (for root process)
  int *students_grade, *city_highest, *city_lowest, *region_highest, *region_lowest;
  float *city_median, *region_median;
  double *city_mean, *city_standard_deviation, *region_mean, *region_standard_deviation;

  float country_median;
  double country_mean, country_standard_deviation;
  int country_highest, country_lowest;

  // variables containing a subset of data (for worker process)
  int *students_grade_subset, *city_highest_subset, *city_lowest_subset, *region_highest_subset, *region_lowest_subset;
  float *city_median_subset, *region_median_subset;
  double *city_mean_subset, *city_standard_deviation_subset, *region_mean_subset, *region_standard_deviation_subset;


  // init input in root process
  if (rank == 0) {
    scanf("%d %d %d %d", &size_regions, &size_cities, &size_students, &seed);
    // allocate root memmory
    students_grade = malloc(sizeof(int) * size_regions * size_cities * size_students);

    city_median = malloc(sizeof(float *) * size_regions * size_cities);
    city_mean = malloc(sizeof(double *) * size_regions * size_cities);
    city_standard_deviation = malloc(sizeof(double *) * size_regions * size_cities);
    city_highest = malloc(sizeof(int *) * size_regions * size_cities);
    city_lowest = malloc(sizeof(int *) * size_regions * size_cities);

    region_median = (float *) malloc(sizeof(float) * size_regions);
    region_mean = (double *) malloc(sizeof(double) * size_regions);
    region_standard_deviation = (double *) malloc(sizeof(double) * size_regions);
    region_highest = (int *) malloc(sizeof(int) * size_regions);
    region_lowest = (int *) malloc(sizeof(int) * size_regions);

    srand(seed);

    for (int i = 0; i < size_regions; i++) {
      for (int j = 0; j < size_cities; j++) {
        for (int k = 0; k < size_students; k++) {
          students_grade[calculate_index(i, j, k, size_regions, size_cities, size_students)] = rand() % 101;
        }
      }
    }
  }


  MPI_Bcast(&size_cities, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&size_regions, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&size_students, 1, MPI_INT, 0, MPI_COMM_WORLD);

  int *displs = malloc(sizeof(int) * comm_size);
  int *send_counts = malloc(sizeof(int) * comm_size);
  int *size_subset = malloc(sizeof(int) * comm_size);
  calculate_counts_displs(send_counts, displs, size_subset, comm_size, size_regions,
                          size_cities * size_students);

  chunk_size = size_subset[rank];
  response_time = omp_get_wtime();
  printf("%d - %d \n", rank, chunk_size);

  students_grade_subset = malloc(sizeof(int) * chunk_size * size_cities * size_students);

  city_median_subset = malloc(sizeof(float *) * chunk_size * size_cities);
  city_mean_subset = malloc(sizeof(double *) * chunk_size * size_cities);
  city_standard_deviation_subset = malloc(sizeof(double *) * chunk_size * size_cities);
  city_highest_subset = malloc(sizeof(int *) * chunk_size * size_cities);
  city_lowest_subset = malloc(sizeof(int *) * chunk_size * size_cities);

  region_median_subset = (float *) malloc(sizeof(float) * chunk_size);
  region_mean_subset = (double *) malloc(sizeof(double) * chunk_size);
  region_standard_deviation_subset = (double *) malloc(sizeof(double) * chunk_size);
  region_highest_subset = (int *) malloc(sizeof(int) * chunk_size);
  region_lowest_subset = (int *) malloc(sizeof(int) * chunk_size);

  // print calculated send counts and displacements for each process
  if (0 == rank) {
    for (int i = 0; i < comm_size; i++) {

    }
  }


  // scatter values to each process
  MPI_Scatterv(students_grade, send_counts, displs, MPI_INT, students_grade_subset,
               send_counts[rank], MPI_INT, 0, MPI_COMM_WORLD);

  if (rank == 0) {
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

  } else {
    calculate_city_stats(
      chunk_size,
      size_cities,
      size_students,
      students_grade_subset,
      city_median_subset,
      city_mean_subset,
      city_highest_subset,
      city_lowest_subset,
      city_standard_deviation_subset
    );

    calculate_region_stats(
      chunk_size,
      size_cities,
      size_students,
      students_grade_subset,
      region_median_subset,
      region_mean_subset,
      region_highest_subset,
      region_lowest_subset,
      region_standard_deviation_subset
    );
  }

  calculate_counts_displs(send_counts, displs, size_subset, comm_size, size_regions,
                          chunk_size * size_students);

//
//  MPI_Gatherv(
//    city_median_subset,
//    chunk_size * size_cities,
//    MPI_FLOAT,
//    city_median,
//    chunk_size * size_cities,
//    MPI_FLOAT,
//    0,
//    MPI_COMM_WORLD);
//
//  MPI_Gather(
//    city_mean_subset,
//    chunk_size * size_cities,
//    MPI_DOUBLE,
//    city_mean,
//    chunk_size * size_cities,
//    MPI_DOUBLE,
//    0,
//    MPI_COMM_WORLD);
//
//  MPI_Gather(
//    city_highest_subset,
//    chunk_size * size_cities,
//    MPI_INT,
//    city_highest,
//    chunk_size * size_cities,
//    MPI_INT,
//    0,
//    MPI_COMM_WORLD);
//
//  MPI_Gather(
//    city_lowest_subset,
//    chunk_size * size_cities,
//    MPI_INT,
//    city_lowest,
//    chunk_size * size_cities,
//    MPI_INT,
//    0,
//    MPI_COMM_WORLD);
//
//  MPI_Gather(
//    city_standard_deviation_subset,
//    chunk_size * size_cities,
//    MPI_DOUBLE,
//    city_standard_deviation,
//    chunk_size * size_cities,
//    MPI_DOUBLE,
//    0,
//    MPI_COMM_WORLD);
//
//  MPI_Gather(
//    region_median_subset,
//    chunk_size,
//    MPI_FLOAT,
//    region_median,
//    chunk_size,
//    MPI_FLOAT,
//    0,
//    MPI_COMM_WORLD);
//
//  MPI_Gather(
//    region_mean_subset,
//    chunk_size,
//    MPI_DOUBLE,
//    region_mean,
//    chunk_size,
//    MPI_DOUBLE,
//    0,
//    MPI_COMM_WORLD);
//
//  MPI_Gather(
//    region_highest_subset,
//    chunk_size,
//    MPI_INT,
//    region_highest,
//    chunk_size,
//    MPI_INT,
//    0,
//    MPI_COMM_WORLD);
//
//  MPI_Gather(
//    region_lowest_subset,
//    chunk_size,
//    MPI_INT,
//    region_lowest,
//    chunk_size,
//    MPI_INT,
//    0,
//    MPI_COMM_WORLD);
//
//  MPI_Gather(
//    region_standard_deviation_subset,
//    chunk_size,
//    MPI_DOUBLE,
//    region_standard_deviation,
//    chunk_size,
//    MPI_DOUBLE,
//    0,
//    MPI_COMM_WORLD);

  // calculate country status in root process because we need the entire data set to be calculated at same time
  if (rank == 0) {


    response_time = omp_get_wtime() - response_time;
    printf("%f", response_time);

    // print results
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

//    free(students_grade);
//
//    free(city_median);
//    free(city_mean);
//    free(city_standard_deviation);
//    free(city_highest);
//    free(city_lowest);
//
//    free(region_median);
//    free(region_mean);
//    free(region_standard_deviation);
//    free(region_highest);
//    free(region_lowest);
  }

  free(students_grade_subset);

  free(city_median_subset);
  free(city_mean_subset);
  free(city_standard_deviation_subset);
  free(city_highest_subset);
  free(city_lowest_subset);

  free(region_median_subset);
  free(region_mean_subset);
  free(region_standard_deviation_subset);
  free(region_highest_subset);
  free(region_lowest_subset);

  MPI_Finalize();

  return 0;
}