#pragma once
#include <string>

int i4_choose ( int n, int k );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int i4_uniform_ab ( int a, int b, int &seed );
void i4vec_print ( int n, int a[], std::string title );
int i4vec_sum ( int n, int a[] );
int *i4vec_uniform_ab_new ( int n, int a, int b, int &seed );
int mono_between_enum ( int m, int n1, int n2 );
void mono_between_next_grevlex ( int m, int n1, int n2, int x[] );
void mono_between_next_grlex ( int m, int n1, int n2, int x[] );
int *mono_between_random ( int m, int n1, int n2, int &seed, int &rank );
void mono_next_grevlex ( int m, int x[] );
void mono_next_grlex ( int m, int x[] );
void mono_print ( int m, int f[], std::string title );
int mono_rank_grlex ( int m, int x[] );
int mono_total_enum ( int m, int n );
void mono_total_next_grevlex ( int m, int n, int x[] );
void mono_total_next_grlex ( int m, int n, int x[] );
int *mono_total_random ( int m, int n, int &seed, int &rank );
int *mono_unrank_grlex ( int m, int rank );
int mono_upto_enum ( int m, int n );
void mono_upto_next_grevlex ( int m, int n, int x[] );
void mono_upto_next_grlex ( int m, int n, int x[] );
int *mono_upto_random ( int m, int n, int &seed, int &rank );
double *mono_value ( int m, int n, int f[], double x[] );
void timestamp ( );
