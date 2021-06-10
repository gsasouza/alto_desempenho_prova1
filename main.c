#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// compare function for qsort
int Compare(const void *a, const void *b) {
  int int_a = *((short int *) a);
  int int_b = *((short int *) b);

  if (int_a == int_b) return 0;
  else if (int_a < int_b) return -1;
  else return 1;
}

short int calculate_lowest(short int *linha, int tam) {
  short int menor = 101;
  for (int i = 0; i < tam; i++) {
    if (linha[i] < menor) {
      menor = linha[i];
    }
  }
  return menor;
}

short int calculate_highest(short int *linha, int tam) {
  short int maior = -1;
  for (int i = 0; i < tam; i++) {
    if (linha[i] > maior) {
      maior = linha[i];
    }
  }
  return maior;
}

float calculate_median(short int *linha, int tam) {
  float mediana;

  qsort(linha, tam, sizeof(int), Compare);

  if (tam % 2) {
    mediana = linha[(int) (tam / 2)];
  } else {
    mediana = (linha[tam / 2 - 1] + linha[tam / 2]) / 2;
  }
  return mediana;
}

float calculate_mean(short int *linha, int tam) {
  int total = 0;
  for (int i = 0; i < tam; i++) {
    total += linha[i];
  }
  return (float) total / (float) tam;
}

// standard deviation
float calculate_standard_deviation(short int *linha, int tam) {
  float media, dp = 0.0;
  int i;

  media = calculate_mean(linha, tam);

  for (i = 0; i < tam; ++i) {
    dp += pow(linha[i] - media, 2);
  }
  return sqrt(dp / (float) tam);
}

void calculate_city_stats(int size_regions, int size_students, int size_cities, short int ***students_grade, float **city_median,
                          float **city_mean,
                          short int **city_highest, short int **city_lowest, float **city_standard_deviation) {
  for (int i = 0; i < size_regions; i++) {
    for (int j = 0; j < size_cities; j++) {
      city_median[i][j] = calculate_median(students_grade[i][j], size_students);
      city_mean[i][j] = calculate_mean(students_grade[i][j], size_students);
      city_highest[i][j] = calculate_highest(students_grade[i][j], size_students);
      city_lowest[i][j] = calculate_lowest(students_grade[i][j], size_students);
      city_standard_deviation[i][j] = calculate_standard_deviation(students_grade[i][j], size_students);
    }
  }
}

short int ***allocate_3d_array(int size_x, int size_y, int size_z) {
  short int ***array = (short int ***) malloc(sizeof(short int ***) * (size_x + 1) );
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
      printf("%f.2 ", matrix[i][j]);
    }
    printf("\n");
  }
}


int main() {
  int seed, size_cities, size_regions, size_students;
  scanf("%d %d %d %d", &size_regions, &size_students, &size_cities, &seed);
  short int ***students_grade = allocate_3d_array(size_regions, size_students, size_cities);
  float **city_median = allocate_2d_array_float(size_regions, size_students);
  float **city_mean = allocate_2d_array_float(size_regions, size_students);
  short int **city_highest = allocate_2d_array(size_regions, size_students);
  short int **city_lowest = allocate_2d_array(size_regions, size_students);
  float **city_standard_deviation = allocate_2d_array_float(size_regions, size_students);
  float *region_median = (float *) malloc(sizeof(float) * size_regions);
  float *region_mean = (float *) malloc(sizeof(float) * size_regions);
  short int *region_highest = (short int *) malloc(sizeof(short int) * size_regions);
  short int *region_lowest = (short int *) malloc(sizeof(short int) * size_regions);
  float *region_standard_deviation = (float *) malloc(sizeof(float) * size_regions);
  float country_median, country_mean, country_standard_deviation;
  short int country_highest, country_lowest;
  srand(seed);

  for (int i = 0; i < size_regions; i++) {
    for (int j = 0; j < size_students; j++) {
      for (int k = 0; k < size_cities; k++) {
        students_grade[i][j][k] = abs((short int) rand() % 100);
      }
    }
  }

  for (int i = 0; i < size_regions; i++) {
    printf("\nRegião %d\n", i);
    for (int j = 0; j < size_students; j++) {
      for (int k = 0; k < size_cities; k++) {
        printf("%d ", students_grade[i][j][k]);
      }
      printf("\n");
    }
  }

  calculate_city_stats(size_regions, size_cities, size_students, students_grade, city_median, city_mean, city_highest, city_lowest,
                       city_standard_deviation);

  printf("mediana \n");
  print_float_matrix(city_median, size_regions, size_cities);

  printf("media \n");
  print_float_matrix(city_mean, size_regions,  size_cities);

  printf("desvio \n");
  print_float_matrix(city_standard_deviation, size_regions, size_cities);

  printf("maior \n");
  print_int_matrix(city_highest, size_regions, size_cities);

  printf("menor \n");
  print_int_matrix(city_lowest, size_regions, size_cities);


  return 0;
}
