/**
 * Param implementation.
 * Copyright 2017 MPI-CBG/MPI-PKS Peter Schwede
 */
#include <string.h>
#include <stdio.h>
#include <assert.h>

#include "EmissionTable.h"
#include "Logging.h"
#include "SafeAlloc.h"

#include "Params.h"

#define MEMORYLIMIT 16

bool Params__create(struct Params* self, struct EmissionTable emission_tables[6]) {
  self->dirty = false;
  self->firstexon = false;
  self->lastexon = false;
  self->max_memory = MEMORYLIMIT;

  self->num_stop_codons = 3;
  self->stop_codons = (Literal*) SAFEMALLOC(sizeof(Literal) * 3 * self->num_stop_codons);
  self->stop_codons[0] = LITERAL_T;
  self->stop_codons[1] = LITERAL_G;
  self->stop_codons[2] = LITERAL_A;

  self->stop_codons[3] = LITERAL_T;
  self->stop_codons[4] = LITERAL_A;
  self->stop_codons[5] = LITERAL_A;

  self->stop_codons[6] = LITERAL_T;
  self->stop_codons[7] = LITERAL_A;
  self->stop_codons[8] = LITERAL_G;

  self->split_emissions_acceptor = 99;
  self->split_emissions_donor = 99;

  self->emission_table_4_UNIFORM = &emission_tables[0];
  self->emission_table_16_UNIFORM = &emission_tables[1];
  self->emission_table_64_UNIFORM = &emission_tables[2];
  self->emission_table_64_LAMBDA = &emission_tables[3];
  self->emission_table_61_LAMBDA = &emission_tables[4];
  self->emission_table_61_UNIFORM = &emission_tables[5];

  //strncpy(self->matrices_path_prefix, "CESAR/matrices", PATH_STRING_LENGTH-100);
  strncpy(self->dot, "\0", PATH_STRING_LENGTH-1);
  strncpy(self->clade, "human", PATH_STRING_LENGTH-1);

  //strncpy(self->blosum_file, "blosum_freq62", PATH_STRING_LENGTH-1);

  strncpy(self->eth_file, "extra/tables/%s/eth_codon_sub.txt", PATH_STRING_LENGTH-1);
  strncpy(self->acc_profile, "extra/tables/%s/acc_profile.txt", PATH_STRING_LENGTH-1);
  strncpy(self->first_codon_profile, "extra/tables/%s/firstCodon_profile.txt", PATH_STRING_LENGTH-1);
  strncpy(self->do_profile, "extra/tables/%s/do_profile.txt", PATH_STRING_LENGTH-1);
  strncpy(self->last_codon_profile, "extra/tables/%s/lastCodon_profile.txt", PATH_STRING_LENGTH-1);
  strncpy(self->u12_acc_profile, "extra/tables/%s/u12_acc_profile.txt", PATH_STRING_LENGTH-1);
  strncpy(self->u12_donor_profile, "extra/tables/%s/u12_donor_profile.txt", PATH_STRING_LENGTH-1);
    
  EmissionTable__init(self->emission_table_4_UNIFORM,  1, UNIFORM_DISTRIBUTION);
  EmissionTable__init(self->emission_table_16_UNIFORM, 2, UNIFORM_DISTRIBUTION);
  EmissionTable__init(self->emission_table_64_UNIFORM, 3, UNIFORM_DISTRIBUTION);
  EmissionTable__init(self->emission_table_64_LAMBDA,  3, LAMBDA_DISTRIBUTION);
  EmissionTable__init(self->emission_table_61_LAMBDA,  3, LAMBDA_DISTRIBUTION);
  EmissionTable__init(self->emission_table_61_UNIFORM, 3, UNIFORM_DISTRIBUTION);

  self->no_leading_introns_prob = 0.5;
  self->no_trailing_introns_prob = 0.49;

  // transition probabilities that are currently 1.0:
  self->acc_acc = 0;
  self->acc2_ii = 0;
  self->do_do = 0;
  self->ii1_ii2 = 0;
  self->i1_i2 = 0;
  self->i2_i3 = 0;
  self->c2_c3 = 0;

  self->stop_codon_emission_logodd = Logodd__log(0.00019625991262287552);  // based on 1704/8,682,364

  long double multiple_cd_factors[] = {
    0.432,  // delete 2 codons at once = 0.432 + 1cd_logodd
    0.276,  // delete 3 codons at once ...;
    0.208,  //        4;
    0.164,  //       ...
    0.147,
    0.133,
    0.127,
    0.123,
    0.118   // delete 10 codons at once ...
  };
  self->num_factors = 9;
  memcpy(self->multiple_cd_factors, multiple_cd_factors, sizeof(long double) * self->num_factors);
  long double sum_multiple_cd_factors = 0;
  for (uint8_t i=0; i < self->num_factors; i++) {
    sum_multiple_cd_factors += self->multiple_cd_factors[i];
  }
  
  self->fs_logodd = Logodd__log(0.0001);
  self->ci_logodd = Logodd__log(0.01);
  self->ci2_logodd = Logodd__log(0.2);

  self->skip_acc = self->fs_logodd;
  self->skip_do = self->fs_logodd;

  self->splice_i1 = self->ci_logodd;
  self->splice_nti = self->fs_logodd;
  self->nti_nti_logodd = Logodd__log(0.25);
  self->total_cd_logodd = Logodd__log(0.025);
  
  // factors for splice site
  self->cd_logodd = self->total_cd_logodd - Logodd__log(1.0 + sum_multiple_cd_factors);
  self->cd_acc = Logodd__log(0.012);
  self->cd_do = Logodd__log(0.012);

  self->c3_i1 = self->ci_logodd;
  self->c3_i1_do = self->ci_logodd;
  self->i3_i1 = self->ci2_logodd;
  self->i3_i1_acc = Logodd__log(0.4);
  self->i3_i1_do = Logodd__log(0.4);

  self->c3_js = self->ci2_logodd;

  self->i3_js = Logodd__subexp(0., self->ci2_logodd);

  self->nti_nti = self->nti_nti_logodd;
  self->nti_js = Logodd__subexp(0., self->nti_nti);

  self->e1_e1 = Logodd__log(0.9);  // TODO test more thoroughly.
                           // only matters with long introns

  self->b2_b2 = Logodd__log(0.9);  // TODO test more thoroughly.
                           // only matters with long introns

  self->no_leading_introns_logodd = Logodd__log(self->no_leading_introns_prob);
  self->no_trailing_introns_logodd = Logodd__log(self->no_trailing_introns_prob);

  self->bsd_e2 = self->fs_logodd;
  self->b1_bas = self->fs_logodd;
  self->b2_bas = self->fs_logodd;
  return true;
}


bool Params__set_paths(struct Params* self) {
  char tmp[PATH_STRING_LENGTH];
  sprintf(tmp, self->eth_file,            self->clade);
  strncpy(self->eth_file, tmp, PATH_STRING_LENGTH);
  sprintf(tmp, self->acc_profile,         self->clade);
  strncpy(self->acc_profile, tmp, PATH_STRING_LENGTH);
  sprintf(tmp, self->first_codon_profile, self->clade);
  strncpy(self->first_codon_profile, tmp, PATH_STRING_LENGTH);
  sprintf(tmp, self->do_profile,          self->clade);
  strncpy(self->do_profile, tmp, PATH_STRING_LENGTH);
  sprintf(tmp, self->last_codon_profile,  self->clade);
  strncpy(self->last_codon_profile, tmp, PATH_STRING_LENGTH);
  sprintf(tmp, self->u12_acc_profile,     self->clade);
  strncpy(self->u12_acc_profile, tmp, PATH_STRING_LENGTH);
  sprintf(tmp, self->u12_donor_profile,   self->clade);
  strncpy(self->u12_donor_profile, tmp, PATH_STRING_LENGTH);

  /*
  logv(1, "fasta_file:\t%s",          self->fasta_file);
  logv(1, "eth_file:\t%s",            self->eth_file);
  logv(1, "acc_profile:\t%s",         self->acc_profile);
  logv(1, "first_codon_profile:\t%s", self->first_codon_profile);
  logv(1, "do_profile:\t%s",          self->do_profile);
  logv(1, "last_codon_profile:\t%s",  self->last_codon_profile);
  logv(1, "u12_acc_profile:\t%s",     self->u12_acc_profile);
  logv(1, "u12_donor_profile:\t%s",   self->u12_donor_profile);
  */

  return true;
}


/** Inserts are like matches against random reference sequences.
 *
 * In other words, inserts are treated like a match for any reference and have
 * the same probability for every reference. However each query has a different
 * probability that depends on the substitution matrix.
 * 
 * We use the substitution matrix to calculate the insert matrix. For that we
 * need two prerequisites:
 *   1. the sum of all probabilities *scores_sum* to normalize the values.
 *   2. the sum of all stop codon emission probabilites *fix_factor* to
 *      re-adjust all non-stop codon probs.
 * Afterwards (3.) we will remove the stop codons by setting all their
 * emission probabilites to zero while we slightly increase all other values
 * by the factor (1+fix_factor) to make sure, those sum up to one again.
 *
 **/
bool Params__make_insert_table(struct Params* self, struct EmissionTable* table) {

  // 1. get the total sum of all probabilities
  LOGODD_T scores_sum = LOGODD_NEGINF, stop_codon_sum = LOGODD_NEGINF;
  for (uint8_t query=0; query < table->values->num_columns; query++) {
    for (uint8_t reference=0; reference < table->values->num_rows; reference++) {
      scores_sum = Logodd__sumexp(scores_sum, EmissionTable__get(table, query, reference));
    }
  }

  // 2. get the sum of normalized stop codon emission probs
  for (uint8_t i=0; i < self->num_stop_codons; i++) {
    LOGODD_T sum_raw_insertion_score = LOGODD_NEGINF;
    for (uint8_t query=0; query < table->values->num_columns; query++) {
      sum_raw_insertion_score = Logodd__sumexp(sum_raw_insertion_score, EmissionTable__get(table, query, Literal__uint(3, &self->stop_codons[i*3])));
    }
    stop_codon_sum = Logodd__sumexp(stop_codon_sum, sum_raw_insertion_score);
  }
  stop_codon_sum = Logodd__add(stop_codon_sum, -scores_sum);

  // 3. set (query x reference)
  LOGODD_T fix_factor = stop_codon_sum/1.0;
  for (uint8_t query=0; query < table->values->num_columns; query++) {

    // handle stop_codons
    bool is_stop_codon = false;
    for (uint8_t i=0; i < self->num_stop_codons; i++) {
      uint8_t stop_codon_id = Literal__uint(3, &self->stop_codons[i*3]);
      is_stop_codon = query == stop_codon_id;
      if (is_stop_codon) {
        LogoddMatrix__set(table->values, query, stop_codon_id, LOGODD_NEGINF);
        break;
      }
    }
    if (is_stop_codon) {
      continue;
    }

    // handle all other
    LOGODD_T sum_raw_insertion_score = LOGODD_NEGINF;
    for (uint8_t reference=0; reference < table->values->num_rows; reference++) {
      sum_raw_insertion_score = Logodd__sumexp(sum_raw_insertion_score, EmissionTable__get(table, query, reference));
    }

    double normalized = sum_raw_insertion_score + Logodd__sumexp(0, fix_factor) - scores_sum;

    for (uint8_t reference=0; reference < table->values->num_rows; reference++) {
      LogoddMatrix__set(table->values, query, reference, normalized);
    }
  }
  
  return true;
}


bool Params__recalculate(struct Params* self) {
  long double sum_multiple_cd_factors = 0;
  for (uint8_t i=0; i < self->num_factors; i++) {
    sum_multiple_cd_factors += self->multiple_cd_factors[i];
  }

  EmissionTable__read(self->emission_table_64_LAMBDA, self->eth_file);

  self->cd_logodd = self->total_cd_logodd - Logodd__log(1.0 + sum_multiple_cd_factors);
  self->splice_js = Logodd__subexp(Logodd__subexp(0, self->fs_logodd), self->splice_i1);

  self->js_js = self->cd_logodd;
  self->js_c1 = Logodd__subexp(Logodd__subexp(0, Logodd__sumexp(self->fs_logodd, self->fs_logodd)), self->total_cd_logodd);
  self->bas_sca = Logodd__subexp(0, self->fs_logodd);

  self->nti_nti = self->nti_nti_logodd;
  self->nti_js = Logodd__subexp(0, self->nti_nti);

  self->c3_js = Logodd__subexp(Logodd__subexp(0, self->fs_logodd), self->c3_i1);

  self->i3_i1 = self->ci2_logodd;
  self->i3_js = Logodd__subexp(0, self->ci2_logodd);
  self->i3_js_do = Logodd__subexp(0, self->i3_i1_do);
  self->i3_js_acc = Logodd__subexp(0, self->i3_i1_acc);

  self->b2_acc = Logodd__subexp(Logodd__subexp(0, self->b2_b2), self->skip_acc);
  self->b1_b2 = Logodd__add(Logodd__subexp(0, self->b1_bas), self->no_leading_introns_logodd);
  self->b1_acc = Logodd__add(Logodd__subexp(0, self->b1_bas), self->no_leading_introns_logodd);

  self->js_scd = Logodd__subexp(0, self->fs_logodd);
  self->bsd_do = Logodd__subexp(Logodd__subexp(0, self->skip_do), self->bsd_e2);
  self->bsd_do_id = Logodd__subexp(Logodd__subexp(Logodd__subexp(0, self->skip_do), self->bsd_e2), self->intron_del);
  self->do2_e1 = Logodd__log(1.0 - self->no_trailing_introns_prob);
  self->do2_e2 = self->no_trailing_introns_logodd;

  self->e1_e2 = Logodd__subexp(0, self->e1_e1);

  // assign substitutions to non-/stop codons
  for (uint8_t this=0; this < self->num_stop_codons; this++) {
    uint8_t row = Literal__uint(3, &self->stop_codons[this*3]);

    for (uint8_t column=0; column < self->emission_table_61_UNIFORM->values->num_columns; column++) {
      bool overwrite = true;
      for (uint8_t other=0; overwrite && other < self->num_stop_codons; other++) {
        if (column == Literal__uint(3, &self->stop_codons[other*3])) {
          overwrite = false;
        }
      }

      LOGODD_T logodd = LogoddMatrix__get(self->emission_table_61_UNIFORM->values, column, row);
      LOGODD_T new_value = self->stop_codon_emission_logodd - Logodd__log(3);
      if (overwrite) {
        new_value = Logodd__subexp(0, Logodd__exp(self->stop_codon_emission_logodd)) + logodd;
      }
      logv(7, "(%ix%i)\t%f -> %f", column, row, Logodd__exp(logodd), Logodd__exp(new_value));
      LogoddMatrix__set(self->emission_table_61_UNIFORM->values, column, row, new_value);
    }
  }

  EmissionTable__read(self->emission_table_61_LAMBDA, self->eth_file);
  Params__make_insert_table(self, self->emission_table_61_LAMBDA);

  if (g_loglevel > 3) {
    FILE* matrix_log = fopen("cesar_matrix.log", "a");
    char tmp[1024000] = "";
    EmissionTable__str(self->emission_table_64_LAMBDA, tmp);
    fprintf(matrix_log, "64Lambda:\n%s", tmp);
    EmissionTable__str(self->emission_table_61_LAMBDA, tmp);
    fprintf(matrix_log, "61Lambda:\n%s", tmp);
    EmissionTable__str(self->emission_table_64_UNIFORM, tmp);
    fprintf(matrix_log, "64Uniform:\n%s", tmp);
    EmissionTable__str(self->emission_table_61_UNIFORM, tmp);
    fprintf(matrix_log, "61Uniform:\n%s", tmp);
    fclose(matrix_log);
  }

  return true;
}

bool Params__set_via_str(struct Params* self, char* string, char* value) {
  #define STRING_LENGTH 10
  const void* STRING_DICT[STRING_LENGTH][2] = {
    {"clade",                     &self->clade},
    {"eth_file",                  &self->eth_file},
    {"acc_profile",               &self->acc_profile},
    {"do_profile",                &self->do_profile},
    {"first_codon_profile",       &self->first_codon_profile},
    {"last_codon_profile",        &self->last_codon_profile},
    {"u12_acc_profile",           &self->u12_acc_profile},
    {"u12_donor_profile",         &self->u12_donor_profile},
    {"fasta_file",                &self->fasta_file},
    {"dot",                       &self->dot}
  };
  #define LODD_LENGTH 51
  const void* LODD_DICT[LODD_LENGTH][2] = {
    {"fs_prob",                   &self->fs_logodd},
    {"ci_prob",                   &self->ci_logodd},
    {"ci2_prob",                  &self->ci2_logodd},
    {"nti_nti_prob",              &self->nti_nti_logodd},
    {"total_cd_prob",             &self->total_cd_logodd},
    {"cd_acc",                    &self->cd_acc},
    {"cd_do",                     &self->cd_do},
    {"c3_i1_do",                  &self->c3_i1_do},
    {"i3_i1_do",                  &self->i3_i1_do},
    {"i3_i1_acc",                 &self->i3_i1_acc},
    {"i3_js_do",                  &self->i3_js_do},
    {"i3_js_acc",                 &self->i3_js_acc},
    {"stop_codon_emission_prob",  &self->stop_codon_emission_logodd},
    {"no_leading_introns_prob",   &self->no_leading_introns_logodd},
    {"no_traling_introns_prob",   &self->no_trailing_introns_logodd},
    {"intron_del",                &self->intron_del},
    {"split_do1",                 &self->split_do1},
    {"js_js",                     &self->js_js},
    {"bas_sca",                   &self->bas_sca},
    {"bsd_do",                    &self->bsd_do},
    {"b1_bas",                    &self->b1_bas},
    {"b2_bas",                    &self->b2_bas},
    {"b1_acc",                    &self->b1_acc},
    {"b1_b2",                     &self->b1_b2},
    {"b2_b2",                     &self->b2_b2},
    {"b2_acc",                    &self->b2_acc},
    {"acc_acc",                   &self->acc_acc},
    {"acc2_ii",                   &self->acc2_ii},
    {"skip_acc",                  &self->skip_acc},
    {"js_scd",                    &self->js_scd},
    {"do_do",                     &self->do_do},
    {"do2_e1",                    &self->do2_e1},
    {"do2_e2",                    &self->do2_e2},
    {"skip_do",                   &self->skip_do},
    {"splice_js",                 &self->splice_js},
    {"splice_nti",                &self->splice_nti},
    {"nti_js",                    &self->nti_js},
    {"nti_nti",                   &self->nti_nti},
    {"splice_i1",                 &self->splice_i1},
    {"ii1_ii2",                   &self->ii1_ii2},
    {"i1_i2",                     &self->i1_i2},
    {"i2_i3",                     &self->i2_i3},
    {"i3_i1",                     &self->i3_i1},
    {"i3_js",                     &self->i3_js},
    {"js_c1",                     &self->js_c1},
    {"c2_c3",                     &self->c2_c3},
    {"c3_i1",                     &self->c3_i1},
    {"c3_js",                     &self->c3_js},
    {"bsd_e2",                    &self->bsd_e2},
    {"e1_e1",                     &self->e1_e1},
    {"e1_e2",                     &self->e1_e2},
  };
#define UINT_LENGTH 3
  const void* UINT_DICT[UINT_LENGTH][2] = {
    {"split_emissions_acceptor",  &self->split_emissions_acceptor},
    {"split_emissions_donor",     &self->split_emissions_donor},
    {"max_memory",                &self->max_memory}
  };

  // assume value is string
  for (int i=0; i < STRING_LENGTH; i++) {
    if (!strcmp(STRING_DICT[i][0], string)) {
      strncpy((char*) STRING_DICT[i][1], value, PATH_STRING_LENGTH);
      logv(1, "Setting %s := %s", string, (char*) STRING_DICT[i][1]);
      return true;
    }
  }

  // assume value is a float
  for (int i=0; i < LODD_LENGTH; i++) {
    if (!strcmp(LODD_DICT[i][0], string)) {
      double prob;
      sscanf(value, "%lf", &prob);
      *((LOGODD_T*) LODD_DICT[i][1]) = Logodd__log(prob);
      logv(1, "Setting %s := %E", string, *((LOGODD_T*) LODD_DICT[i][1]));
      self->dirty = true;
      return true;
    }
  }

  // assume value is a uint
  for (int i=0; i < UINT_LENGTH; i++) {
    if (!strcmp(UINT_DICT[i][0], string)) {
      unsigned int uint = 0;
      sscanf(value, "%u", &uint);
      *((uint8_t*) UINT_DICT[i][1]) = (uint8_t) uint;
      logv(1, "Setting %s := %u", string, *((uint8_t*) UINT_DICT[i][1]));
      return true;
    }
  }

  return false;
}

void Params__destroy(struct Params* self) {
  EmissionTable__destroy(self->emission_table_4_UNIFORM);
  EmissionTable__destroy(self->emission_table_16_UNIFORM);
  EmissionTable__destroy(self->emission_table_64_UNIFORM);
  EmissionTable__destroy(self->emission_table_64_LAMBDA);
  EmissionTable__destroy(self->emission_table_61_LAMBDA);
  EmissionTable__destroy(self->emission_table_61_UNIFORM);
  free(self->stop_codons);
}