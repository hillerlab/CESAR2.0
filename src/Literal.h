/**
 * Literal definition.
 * Copyright 2017 MPI-CBG/MPI-PKS Peter Schwede
 */
#ifndef LITERAL_H_
#define LITERAL_H_

#include <stddef.h>
#include <stdint.h>

#define EMISSION_ID_T uint_fast8_t

typedef enum Literal {
  LITERAL_A,
  LITERAL_C,
  LITERAL_G,
  LITERAL_T,
  LITERAL_N
} Literal;
#define NUM_LITERALS 4

Literal Literal__from_char(char c);
char Literal__char(Literal literal);
void Literal__str(size_t length, Literal literals[length], char buffer[length]);
EMISSION_ID_T Literal__uint(uint8_t length, Literal array[length]);
uint8_t Literal__Ns(uint8_t length, Literal array[length]);

#endif  // LITERAL_H_
