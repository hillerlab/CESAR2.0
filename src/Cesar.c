/**
 * CESAR entry point.
 * Copyright 2017 MPI-CBG/MPI-PKS Peter Schwede
 */

#include <stdio.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include "Fasta.h"
#include "Profile.h"
#include "Logging.h"
#include "SafeAlloc.h"
#include "Viterbi.h"
#include "Alignment.h"
#include "Arguments.h"
#include "EmissionTable.h"
#include "Sequence.h"
#include "Model.h"


char g_loglevel = LOGLEVEL;

/**
 * The main program.
 */
int main(int argc, char* argv[argc]) {
  struct EmissionTable emission_tables[8];

  struct Params parameters;
  Params__create(&parameters, emission_tables);

  if (!Arguments__read(argc, argv, &parameters)) {
    return 1;
  }

  Params__set_paths(&parameters);

  struct Fasta fasta;
  Fasta__init(&fasta);
  Fasta__read(&fasta, parameters.fasta_file);

  struct Profile** acceptors = SAFEMALLOC(sizeof(Profile*) * fasta.num_references);
  struct Profile** donors = SAFEMALLOC(sizeof(Profile*) * fasta.num_references);

  char prefix[PATH_STRING_LENGTH] = "extra/tables/";

  for (uint8_t i=0; i < fasta.num_references; i++) {
    struct Sequence* reference = fasta.references[i];

    if (reference->acceptor[0] != '\0') {
      char name[STATE_NAME_LENGTH];
      sprintf(name, "ref%iacc", i);
      acceptors[i] = Profile__create(name);
      char path[PROFILE_FILENAME_LENGTH];
      sprintf(path, "%s%s/%s", prefix, parameters.clade, reference->acceptor);
      if (access(path, R_OK) != -1) {
        Profile__read(acceptors[i], path);
      } else {
        Profile__read(acceptors[i], reference->acceptor);
      }
    } else {
      if (fasta.num_references > 1) {
        warn("Missing acceptor profile for reference %u.", i);
      }
      acceptors[i] = Profile__create("acceptor");
		/* only read the default extra/tables/{clade}/firstCodon_profile.txt if -p {acc_profile} {do_profile} is not given (flag parameters.accSpecified) */
      if (i == 0 && !parameters.acc_do_specified && ((parameters.firstexon && fasta.num_references == 1) || (!parameters.firstexon && fasta.num_references > 1))) {
        Profile__read(acceptors[i], parameters.first_codon_profile);
      } else {
        Profile__read(acceptors[i], parameters.acc_profile);
      }
    }
    logv(1, "Reference %u uses acceptor:\t%s", i, acceptors[i]->filename);

    if (reference->acceptor[0] != '\0') {
      char name[STATE_NAME_LENGTH];
      sprintf(name, "ref%idon", i);
      donors[i] = Profile__create(name);
      char path[PROFILE_FILENAME_LENGTH];
      sprintf(path, "%s%s/%s", prefix, parameters.clade, reference->donor);
      if (access(path, R_OK) != -1) {
        Profile__read(donors[i], path);
      } else {
        Profile__read(donors[i], reference->donor);
      }
    } else {
      if (fasta.num_references > 1) {
        warn("Missing donor profile for reference %u.", i);
      }
      donors[i] = Profile__create("donor");
		/* only read the default extra/tables/{clade}/firstCodon_profile.txt if -p {acc_profile} {do_profile} is not given (flag parameters.doSpecified) */
      if (i == fasta.num_references-1 && !parameters.acc_do_specified && ((parameters.lastexon && fasta.num_references == 1) || (!parameters.lastexon && fasta.num_references > 1))) {
        Profile__read(donors[i], parameters.last_codon_profile);
      } else {
        Profile__read(donors[i], parameters.do_profile);
      }
    }
    logv(1, "Reference %u uses donor:\t%s", i, donors[i]->filename);
  }


  /* estimate memory consumption */
  size_t num_states = 0;
  size_t rlength = 0;
  size_t qlength = 0;
  size_t qlength_max = 0;
  for (uint8_t i=0; i < fasta.num_references; i++) {
    struct Sequence* reference = fasta.references[i];
    num_states += 6 + 6 * reference->num_codons + 1 + 2 + 2 + 22 + 6;  /* 22 and 6 for acc and donor states */

    logv(1, "Reference %u length: %lu", i, fasta.references[i]->length);
    logv(1, "Reference %u split codon lengths: %u %u", i, fasta.references[i]->start_split_length, fasta.references[i]->end_split_length);
    rlength += fasta.references[i]->length;
  }
  for (uint8_t i=0; i < fasta.num_queries; i++) {
    logv(1, "Query %u length: %lu", i, fasta.queries[i]->length);
    qlength += fasta.queries[i]->length;
    if (fasta.queries[i]->length > qlength_max) {
	    qlength_max = fasta.queries[i]->length;
    }
  }

  double mem = 
     (num_states * 4 * sizeof(double))                     +   /* Viterbi logodds vmatrix, double has typically 8 bytes */
     (num_states * qlength_max * sizeof(uint32_t))         +   /* Viterbi path matrix, uint32 = 4 bytes*/
     (num_states * sizeof(State))                          +   /* Space for all the state structures, each 304 bytes for a struct State */
	  ((2 * qlength_max + rlength) * sizeof(struct State*)) +   /* Space for Viterbi path, 8 bytes for a pointer on a 64 bit system */
	  ((qlength_max + rlength) * 2 * sizeof(char))          +   /* Space for the aligned seqs */
	  ((qlength_max + rlength) * 1)                         +   /* Space to store all the sequences */
	  1000000;                                                  /* need some smaller emission tables (32 k) and other small stuff */
  mem = mem/1000000000;  /* in GB */	  
  logv(1, "Expecting a memory consumption of: %f GB (max_memory %f)", mem, (double)parameters.max_memory);
  if (mem > (double)parameters.max_memory) {
    die("The memory consumption is limited to %1.4f GB by default. Your attempt requires %1.4f GB. You can change the limit via --max-memory.", (double)parameters.max_memory, mem);
  } else {
    logv(1, "Expecting a memory consumption of: %1.4f GB", mem);
  }


  if (g_loglevel >= 7) {
    char* tmp;
    for (uint8_t i=0; i < fasta.num_references; i++) {
      tmp = SAFECALLOC(sizeof(char), SEQUENCENAMELENGTH + fasta.references[i]->length);
      Literal__str(fasta.references[i]->length, fasta.references[i]->sequence, tmp);
      logv(1, ">original %s\n%s", fasta.references[i]->name, tmp);
      free(tmp);
    }

    for (uint8_t i=0; i < fasta.num_queries; i++) {
      tmp = SAFECALLOC(sizeof(char), SEQUENCENAMELENGTH + fasta.queries[i]->length);
      Literal__str(fasta.queries[i]->length, fasta.queries[i]->sequence, tmp);
      logv(1, ">original %s\n%s", fasta.queries[i]->name, tmp);
      free(tmp);
    }
  }

  clock_t time = clock();
  struct HMM* hmm = multi_exon(&parameters, &fasta, acceptors, donors);
  logv(1, "HMM construction (sec):\t%f", (float)(clock() - time) / CLOCKS_PER_SEC);

  if (parameters.dot[0] != '\0') {
    FILE * dotfile = fopen(parameters.dot, "w");
    HMM__dot(hmm, dotfile);
    fclose(dotfile);
  }

  for (uint8_t q=0; q < fasta.num_queries; q++) {

    size_t length = 2*fasta.queries[q]->length;
    for (uint8_t i=0; i < fasta.num_references; i++) {
      length += fasta.references[i]->length;
    }

    struct State** path = (struct State**) SAFEMALLOC(sizeof(struct State*) * length);
    size_t path_length = 0;

    time = clock();
    Viterbi(hmm, fasta.queries[q]->length, fasta.queries[q]->sequence, &path_length, path);
    logv(1, "Viterbi (sec):\t%f", (float)(clock() - time) / CLOCKS_PER_SEC);

    struct Alignment* alignment = Alignment__create(&fasta, q, &parameters, path_length, path);
    printf(">%s\n%s\n", "referenceExon", alignment->reference);
    printf(">%s\n%s\n", fasta.queries[q]->name, alignment->query);
    Alignment__destroy(alignment);

    free(path);

  }

  HMM__destroy(hmm);
  for (uint8_t i=0; i<fasta.num_references; i++) {
    if(acceptors[i] != NULL) {
      Profile__destroy(acceptors[i]);
    }
    if(donors[i] != NULL) {
      Profile__destroy(donors[i]);
    }
  }
  free(acceptors);
  free(donors);
  Fasta__destroy(&fasta);
  Params__destroy(&parameters);
  return 0;
}
