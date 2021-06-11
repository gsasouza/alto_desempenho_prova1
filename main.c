#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int compare_short(const void *a, const void *b) {
  short int int_a = *((short int *) a);
  short int int_b = *((short int *) b);

  if (int_a == int_b) return 0;
  else if (int_a < int_b) return -1;
  else return 1;
}

short int calculate_lowest_short(short int *linha, int tam) {
  short int menor = 101;
  for (int i = 0; i < tam; i++) {
    if (linha[i] < menor) {
      menor = linha[i];
    }
  }
  return menor;
}

short int calculate_highest_short(short int *linha, int tam) {
  short int maior = -1;
  for (int i = 0; i < tam; i++) {
    if (linha[i] > maior) {
      maior = linha[i];
    }
  }
  return maior;
}

float calculate_median_short(short int *linha, int tam) {
  float mediana;

  qsort(linha, tam, sizeof(short int), compare_short);

  if (tam % 2) {
    mediana = linha[(int) (tam / 2)];
  } else {
    mediana = (linha[tam / 2 - 1] + linha[tam / 2]) / 2;
  }
  return mediana;
}

float calculate_mean_short(short int *linha, int tam) {
  int total = 0;
  for (int i = 0; i < tam; i++) {
    total += linha[i];
  }
  return (float) total / (float) tam;
}

// standard deviation
float calculate_standard_deviation_short(short int *linha, int tam) {
  float media, dp = 0.0;
  int i;

  media = calculate_mean_short(linha, tam);

  for (i = 0; i < tam; ++i) {
    dp += pow(linha[i] - media, 2);
  }
  return sqrt(dp / (float) tam);
}

float calculate_median_float(float *linha, int tam) {
  float mediana;

  qsort(linha, tam, sizeof(short int), compare_short);

  if (tam % 2) {
    mediana = linha[(int) (tam / 2)];
  } else {
    mediana = (linha[tam / 2 - 1] + linha[tam / 2]) / 2;
  }
  return mediana;
}

int compare_float(const void *a, const void *b) {
  float int_a = *((float *) a);
  float int_b = *((float *) b);

  if (int_a == int_b) return 0;
  else if (int_a < int_b) return -1;
  else return 1;
}

float calculate_lowest_float(float *linha, int tam) {
  float menor = 101;
  for (int i = 0; i < tam; i++) {
    if (linha[i] < menor) {
      menor = linha[i];
    }
  }
  return menor;
}

float calculate_mean_float(float *linha, int tam) {
  float total = 0;
  for (int i = 0; i < tam; i++) {
    total += linha[i];
  }
  return (float) total / (float) tam;
}

// standard deviation
float calculate_standard_deviation_float(float *linha, int tam) {
  float media, dp = 0.0;
  int i;

  media = calculate_mean_float(linha, tam);

  for (i = 0; i < tam; ++i) {
    dp += pow(linha[i] - media, 2);
  }
  return sqrt(dp / (float) tam);
}

void calculate_city_stats(int size_regions, int size_cities, int size_students, short int ***students_grade,
                          float **city_median,
                          float **city_mean,
                          short int **city_highest, short int **city_lowest, float **city_standard_deviation) {

  for (int i = 0; i < size_regions; i++) {
    for (int j = 0; j < size_cities; j++) {
      city_median[i][j] = calculate_median_short(students_grade[i][j], size_students);
      city_mean[i][j] = calculate_mean_short(students_grade[i][j], size_students);
      city_highest[i][j] = calculate_highest_short(students_grade[i][j], size_students);
      city_lowest[i][j] = calculate_lowest_short(students_grade[i][j], size_students);
      city_standard_deviation[i][j] = calculate_standard_deviation_short(students_grade[i][j], size_students);
    }
  }
}

void calculate_region_stats(
  int size_regions,
  int size_cities,
  float **city_median,
  float **city_mean,
  short int **city_highest,
  short int **city_lowest,
  float **city_standard_deviation,
  float *region_median,
  float *region_mean,
  short int *region_highest,
  short int *region_lowest,
  float *region_standard_deviation
) {

  for (int i = 0; i < size_regions; i++) {
    region_median[i] = calculate_median_float(city_median[i], size_cities);
    region_mean[i] = calculate_mean_float(city_mean[i], size_cities);
    region_highest[i] = calculate_highest_short(city_highest[i], size_cities);
    region_lowest[i] = calculate_lowest_short(city_lowest[i], size_cities);
    region_standard_deviation[i] = calculate_standard_deviation_float(city_standard_deviation[i], size_cities);
  }
}

void calculate_country_stats(
  int size_regions,
  float *region_median,
  float *region_mean,
  short int *region_highest,
  short int *region_lowest,
  float *region_standard_deviation,
  float *country_median,
  float *country_mean,
  short int *country_highest,
  short int *country_lowest,
  float *country_standard_deviation
) {
  *country_median = calculate_median_float(region_median, size_regions);
  *country_mean = calculate_mean_float(region_mean, size_regions);
  *country_highest = calculate_highest_short(region_highest, size_regions);
  *country_lowest = calculate_lowest_short(region_lowest, size_regions);
  *country_standard_deviation = calculate_standard_deviation_float(region_standard_deviation, size_regions);
}

short int ***allocate_3d_array(int size_x, int size_y, int size_z) {
  short int ***array = (short int ***) malloc(sizeof(short int ***) * (size_x + 1));
  for (int i = 0; i < size_x; i++) {
    array[i] = (short int **) malloc(sizeof(short int **) * (size_y + 1));
    for (int j = 0; j < size_y; j++) {
      array[i][j] = (short int *) malloc(sizeof(short int) * (size_z + 1));
    }
  }
  return array;
}

short int **allocate_2d_array(int size_x, int size_y) {
  short int **array = (short int **) malloc(sizeof(short int *) * size_x);
  for (int i = 0; i < size_x; i++) {
    array[i] = (short int *) malloc(sizeof(short int) * size_y);

  }
  return array;
}

float **allocate_2d_array_float(int size_x, int size_y) {
  float **array = (float **) malloc(sizeof(float *) * size_x);
  for (int i = 0; i < size_x; i++) {
    array[i] = (float *) malloc(sizeof(float) * size_y);

  }
  return array;
}

void print_int_matrix(short int **matrix, int size_x, int size_y) {
  for (int i = 0; i < size_x; i++) {
    printf("\nRegião %d\n", i);
    for (int j = 0; j < size_y; j++) {
      printf("%d ", matrix[i][j]);
    }
    printf("\n");
  }
}

void print_float_matrix(float **matrix, int size_x, int size_y) {
  for (int i = 0; i < size_x; i++) {
    printf("\nRegião %d\n", i);
    for (int j = 0; j < size_y; j++) {
      printf("%.2f ", matrix[i][j]);
    }
    printf("\n");
  }
}

void pretty_print_output(
  int size_regions,
  int size_cities,
  float **city_median,
  float **city_mean,
  short int **city_highest,
  short int **city_lowest,
  float **city_standard_deviation,
  float country_median,
  float country_mean,
  short int country_highest,
  short int country_lowest,
  float country_standard_deviation
) {
  for (int i = 0; i < size_regions; i++) {
    for (int j = 0; j < size_cities; j++) {
      printf("Reg %d - Cid %d: menor: %d, maior: %d, mediana: %.2f, média: %.2f e DP: %.2f\n", i, j, city_lowest[i][j],
             city_highest[i][j], city_median[i][j], city_mean[i][j], city_standard_deviation[i][j]);
    }
    printf("\n\n");
  }
  printf("Brasil: menor: %d, maior: %d, mediana %.2f, média: %.2f e DP: %.2f \n\n", country_lowest, country_highest,
         country_median, country_mean, country_standard_deviation);
}


int main() {
  int seed, size_students, size_regions, size_cities;
  scanf("%d %d %d %d", &size_regions, &size_cities, &size_students, &seed);
  short int ***students_grade = allocate_3d_array(size_regions, size_cities, size_students);
  float **city_median = allocate_2d_array_float(size_regions, size_cities);
  float **city_mean = allocate_2d_array_float(size_regions, size_cities);
  short int **city_highest = allocate_2d_array(size_regions, size_cities);
  short int **city_lowest = allocate_2d_array(size_regions, size_cities);
  float **city_standard_deviation = allocate_2d_array_float(size_regions, size_cities);
  float *region_median = (float *) malloc(sizeof(float) * size_regions);
  float *region_mean = (float *) malloc(sizeof(float) * size_regions);
  short int *region_highest = (short int *) malloc(sizeof(short int) * size_regions);
  short int *region_lowest = (short int *) malloc(sizeof(short int) * size_regions);
  float *region_standard_deviation = (float *) malloc(sizeof(float) * size_regions);
  float country_median, country_mean, country_standard_deviation;
  short int country_highest, country_lowest;
  srand(seed);

  for (int i = 0; i < size_regions; i++) {
    for (int j = 0; j < size_cities; j++) {
      for (int k = 0; k < size_students; k++) {
        students_grade[i][j][k] = rand() % 100;
      }
    }
  }
//
//  for (int i = 0; i < size_regions; i++) {
//    printf("\nRegião %d\n", i);
//    for (int j = 0; j < size_cities; j++) {
//      printf("\n");
//      for (int k = 0; k < size_students; k++) {
//        printf("%d ", students_grade[i][j][k]);
//      }
//    }
//  }

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
    city_median,
    city_mean,
    city_highest,
    city_lowest,
    city_standard_deviation,
    region_median,
    region_mean,
    region_highest,
    region_lowest,
    region_standard_deviation
  );

  calculate_country_stats(
    size_regions,
    region_median,
    region_mean,
    region_highest,
    region_lowest,
    region_standard_deviation,
    &country_median,
    &country_mean,
    &country_highest,
    &country_lowest,
    &country_standard_deviation
  );

//  printf("\nmediana");
//  print_float_matrix(city_median, size_regions, size_cities);
//  printf("\nmedia");
//  print_float_matrix(city_mean, size_regions, size_cities);
//  printf("\nstd");
//  print_float_matrix(city_standard_deviation, size_regions, size_cities);
//  printf("\nmaior");
//  print_int_matrix(city_highest, size_regions, size_cities);
//  printf("\nmenor");
//  print_int_matrix(city_lowest, size_regions, size_cities);
  pretty_print_output(
    size_regions,
    size_cities,
    city_median,
    city_mean,
    city_highest,
    city_lowest,
    city_standard_deviation,
    country_median,
    country_mean,
    country_highest,
    country_lowest,
    country_standard_deviation
  );

  return 0;
}