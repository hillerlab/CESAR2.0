/**
 * Align an exon to a reference
 * Copyright 2017 MPI-CBG/MPI-PKS Peter Schwede
 */

#include <assert.h>
#include <stdbool.h>
#include <string.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "Logging.h"
#include "SafeAlloc.h"
#include "Distribution.h"
#include "EmissionTable.h"
#include "HMM.h"
#include "State.h"
#include "Params.h"
#include "Profile.h"
#include "Fasta.h"

#include "Model.h"


bool create_profile_chain(struct HMM* hmm, struct Profile* profile, struct State** start, struct State** end) {
  struct State* previous_state = NULL;
  for (int i=0; i < profile->length; i++) {
    struct State* current_state = HMM__new_state(hmm);

    Literal reference_base[1];
    reference_base[0] = LITERAL_A;

    char name[STATE_NAME_LENGTH];
    sprintf(name, "match_%s", profile->name);
    State__init(current_state, name, 1, reference_base, &profile->emission_tables[i]);
    current_state->custom = i;
    if (previous_state != NULL) {
      State__add_incoming(previous_state, 0, current_state);
    }
    
    if (i == 0) {
      *start = current_state;
    } else if (i >= profile->length-1) {
      *end = current_state;
    }

    previous_state = current_state;
  }  // Chain end.
  assert(*start != NULL);
  assert(*end != NULL);
  return true;
}


bool forward_deletions(struct HMM* hmm, struct Params* params, size_t num_codons, struct State** codons) {
  for (size_t i=0; i < num_codons; i++) {
    for (uint8_t j=0; j < params->num_factors; j++) {
      if (i+j+1 >= num_codons) {
        break;
      }
      if (i == 0) {
        State__add_incoming(codons[i], params->multiple_cd_factors[j] + params->cd_acc, codons[i+j+2]);
        logv(5, "cd_acc: %e	from: %lu  to: %lu", Logodd__exp(params->multiple_cd_factors[j]+params->cd_acc), i, i+j+2);
      } else if (num_codons == i+j+2) {
        State__add_incoming(codons[i], params->multiple_cd_factors[j] + params->cd_do, codons[i+j+2]);
        logv(5, "cd_do: %e	from: %lu  to: %lu", Logodd__exp(params->multiple_cd_factors[j]+params->cd_do), i, i+j+2);
      } else {
        State__add_incoming(codons[i], Logodd__log(params->multiple_cd_factors[j] * Logodd__exp(params->cd_logodd)), codons[i+j+2]);
      }
    }
#if !NONORMALIZE
    HMM__normalize(hmm, codons[i]);
#endif  // NONORMALIZE
  }
  return true;
}

bool match_codon(struct HMM* hmm, Params* params, size_t index, Literal codon[3], bool last, struct State** curr_cluster) {
  char codon_str[4] = "";
  Literal__str(3, codon, codon_str);

  /*
  for (uint8_t i=0; i < params->num_stop_codons; i++) {
    if (codon[0] == params->stop_codons[i*3] &&
        codon[0+1] == params->stop_codons[i*3+1] &&
        codon[0+2] == params->stop_codons[i*3+2]) {
      die("Reference contains full stop codon %s: Codon %lu", codon_str, index);
    }
  }
  */

  // match_curr_codon
  struct State* match_codon = HMM__new_state(hmm);
  struct State* insert_codon_codon = HMM__new_state(hmm);
  struct State* insert_1nt_codon = HMM__new_state(hmm);
  struct State* delete_1nt_codon = HMM__new_state(hmm);
  struct State* delete_2nt_codon = HMM__new_state(hmm);
  struct State* next_cluster = HMM__new_state(hmm);

  State__init(match_codon, "match_codon", 3, codon, params->emission_table_64_LAMBDA);
  match_codon->custom = index;
  State__init_uniform(insert_codon_codon, "insert_codon_codon", 3, params->emission_table_61_LAMBDA);
  insert_codon_codon->custom = index;
  State__init_uniform(insert_1nt_codon, "insert_1nt_codon", 1, params->emission_table_4_UNIFORM);
  insert_1nt_codon->custom = index;
  State__init_uniform(delete_1nt_codon, "delete_1nt_codon", 2, params->emission_table_16_UNIFORM);
  delete_1nt_codon->custom = index;
  State__init_uniform(delete_2nt_codon, "delete_2nt_codon", 1, params->emission_table_4_UNIFORM);
  delete_2nt_codon->custom = index;
  State__init_silent(next_cluster, "end_codon");
  next_cluster->custom = index;
  
  State__add_incoming(*curr_cluster,      params->js_c1,   match_codon);
  
  State__add_incoming(match_codon,        params->fs_logodd, insert_1nt_codon);
  State__add_incoming(insert_1nt_codon,   params->nti_js,  next_cluster);
  if (last) {
    State__add_incoming(*curr_cluster,      params->cd_do,   next_cluster);
    State__add_incoming(match_codon,        params->c3_i1_do,   insert_codon_codon);
    State__add_incoming(insert_codon_codon, params->i3_i1_do,   insert_codon_codon);
    State__add_incoming(insert_codon_codon, params->i3_js_do,   next_cluster);
  } else {
    if (index == 0) {
      State__add_incoming(*curr_cluster,      params->cd_acc,   next_cluster);
    } else {
      State__add_incoming(*curr_cluster,      params->js_js,   next_cluster);
    }
    State__add_incoming(match_codon,        params->c3_i1,   insert_codon_codon);
    State__add_incoming(insert_codon_codon, params->i3_i1,   insert_codon_codon);
    State__add_incoming(insert_codon_codon, params->i3_js,   next_cluster);
  }
  State__add_incoming(insert_1nt_codon,   params->nti_nti, insert_1nt_codon);

  State__add_incoming(match_codon,        params->c3_js,   next_cluster);



  State__add_incoming(*curr_cluster,      params->fs_logodd, delete_1nt_codon);
  State__add_incoming(*curr_cluster,      params->fs_logodd, delete_2nt_codon);
  State__add_incoming(delete_1nt_codon,   0,               next_cluster);
  State__add_incoming(delete_2nt_codon,   0,               next_cluster);

  *curr_cluster = next_cluster;
  return true;
}


bool create_codon_chain(struct HMM* hmm, struct Params* params, size_t* num_codons, struct State** codons, size_t num_literals, Literal* reference, struct State* start, struct State** end) {
  if (num_literals % 3 != 0) {
    die("A reference is out of frame by %lu nt.", num_literals % 3);
  }
  
  Literal codon[3];
  codons[0] = start;
  *num_codons = 0;

  struct State* curr_cluster = start;
  for (size_t i=0; i+2 < num_literals; i+=3) {
    codon[0] = reference[i];
    codon[1] = reference[i+1];
    codon[2] = reference[i+2];
    //logv(8, "Cluster i=%lu\tnum_codons=%lu\tcodon=%x%x%x", i, (*num_codons)+1, codon[0], codon[1], codon[2]);
    // Create a cluster for this codon.
    match_codon(hmm, params, *num_codons, codon, i+3 >= num_literals, &curr_cluster);
    codons[1 + (*num_codons)++] = curr_cluster;
  }
  *end = curr_cluster;
  return true;
}


struct HMM* multi_exon(struct Params* params, struct Fasta* fasta, struct Profile** acceptors, struct Profile** donors) {
  // read profiles and count states
  size_t num_states = 0;
  for (uint8_t i=0; i < fasta->num_references; i++) {
    struct Sequence* reference = fasta->references[i];

    num_states += 6 + 6 * reference->num_codons + 1 + 2 + 2;

    num_states += acceptors[i]->length;
    num_states += donors[i]->length;
  }

  struct HMM* hmm = HMM__create(num_states, 3, 3);
  Params__recalculate(params);

  uint8_t num_split_codons = 0;
  struct State** split_codons = (struct State**) SAFEMALLOC(sizeof(struct State*) * fasta->num_references * 2);
  struct State* former_intron = NULL;
  struct State* first_intron = NULL;

  logv(1, "There are %i references.", fasta->num_references);
  for (uint8_t i=0; i < fasta->num_references; i++) {
    struct Sequence* reference = fasta->references[i];

    // intron start
    struct State* intron_start;
    if (former_intron == NULL) {
      intron_start = HMM__new_state(hmm);
      State__init_uniform(intron_start, "intron_start", 1, params->emission_table_4_UNIFORM);
      State__add_incoming(intron_start,        params->b2_b2,      intron_start);
      first_intron = intron_start;
    } else {
      intron_start = former_intron;
    }

    // acceptor
    struct State* match_acceptor = NULL;
    struct State* match_acceptor_end = NULL;
    create_profile_chain(hmm, acceptors[i], &match_acceptor, &match_acceptor_end);
    State__add_incoming(intron_start,          params->b2_acc,    match_acceptor);

    struct State* between_acc_split = HMM__new_state(hmm);
    State__init_silent(between_acc_split, "between_acc_split");
    split_codons[num_split_codons++] = between_acc_split;
    State__add_incoming(intron_start,        params->b2_bas,  between_acc_split);
    State__add_incoming(match_acceptor_end,    0,                  between_acc_split);

    logv(1, "reference->start_split_length: %u", reference->start_split_length);
    struct State* split_codon_acceptor = HMM__new_state(hmm);
    switch (reference->start_split_length) {
      case 0:
        State__init_silent(split_codon_acceptor, "split_codon_acceptor");
        break;
      case 1:
        State__init_uniform(split_codon_acceptor, "split_codon_acceptor", 1, params->emission_table_4_UNIFORM);
        break;
      case 2:
        State__init_uniform(split_codon_acceptor, "split_codon_acceptor", 2, params->emission_table_16_UNIFORM);
        break;
      default:
        die("Invalid number of split codon nucleotides: %u", reference->start_split_length);
    }
    State__add_incoming(between_acc_split,     Logodd__log(1.0 - Logodd__exp(params->fs_logodd)),    split_codon_acceptor);

    struct State* insert_1nt_acceptor = HMM__new_state(hmm);
    State__init_uniform(insert_1nt_acceptor, "insert_1nt_acceptor", 1, params->emission_table_4_UNIFORM);
    State__add_incoming(insert_1nt_acceptor,   params->nti_nti,    insert_1nt_acceptor);

    struct State* insert_codon_acceptor = HMM__new_state(hmm);
    State__init_uniform(insert_codon_acceptor, "insert_codon_acceptor", 3, params->emission_table_61_LAMBDA);
    State__add_incoming(insert_codon_acceptor, params->i3_i1_acc,  insert_codon_acceptor);

    struct State* between_split_ins = HMM__new_state(hmm);
    State__init_silent(between_split_ins, "between_split_ins");
    State__add_incoming(split_codon_acceptor,                 0.0,     between_split_ins);
    State__add_incoming(between_acc_split,      params->fs_logodd,     between_split_ins);
    State__add_incoming(between_split_ins,     params->splice_nti,   insert_1nt_acceptor);
    State__add_incoming(between_split_ins,      params->splice_i1, insert_codon_acceptor);

    // codons
    struct State* start_first_codon = HMM__new_state(hmm);
    State__init_silent(start_first_codon, "start_first_codon");
    State__add_incoming(between_split_ins,  params->splice_js,  start_first_codon);
    State__add_incoming(insert_1nt_acceptor,   params->nti_js,     start_first_codon);
    State__add_incoming(insert_codon_acceptor, params->i3_js_acc,      start_first_codon);

    size_t num_codons = 0;
    struct State* end_last_codon = NULL;
    struct State** codons = (struct State**) SAFEMALLOC(sizeof(struct State*) * (2 + reference->num_codons));

    // Iterate over Codons
    create_codon_chain(hmm, params, &num_codons, codons, reference->num_codon_bases, &reference->sequence[reference->codons_offset], start_first_codon, &end_last_codon);

    logv(1, "Codons:\t%lu", num_codons);
    assert(num_codons == reference->num_codons);

    forward_deletions(hmm, params, num_codons, codons);

    free(codons);

    // donor
    logv(1, "reference->end_split_length: %u", reference->end_split_length);
    struct State* split_codon_donor = HMM__new_state(hmm);
    switch (reference->end_split_length) {
      case 0:
        State__init_silent(split_codon_donor, "split_codon_donor");
        break;
      case 1:
        State__init_uniform(split_codon_donor, "split_codon_donor", 1, params->emission_table_4_UNIFORM);
        break;
      case 2:
        State__init_uniform(split_codon_donor, "split_codon_donor", 2, params->emission_table_16_UNIFORM);
        break;
      default:
        die("Invalid number of split codon nucleotides in file %s: %u", params->fasta_file, params->split_emissions_donor);
    }
    State__add_incoming(end_last_codon,        params->js_scd,     split_codon_donor);

    struct State* between_split_donor = HMM__new_state(hmm);
    State__init_silent(between_split_donor, "between_split_donor");
    split_codons[num_split_codons++] = between_split_donor;

    State__add_incoming(end_last_codon,        params->fs_logodd,    between_split_donor);
    State__add_incoming(split_codon_donor,     0.,                   between_split_donor);

    struct State* match_donor = NULL;
    struct State* match_donor_end = NULL;
    create_profile_chain(hmm, donors[i], &match_donor, &match_donor_end);
    if (i+1 == fasta->num_references) {
      State__add_incoming(between_split_donor,   params->bsd_do,      match_donor);
    } else {
      State__add_incoming(between_split_donor,   params->bsd_do_id,   match_donor);
    }

    // intron end
    struct State* intron_end = HMM__new_state(hmm);
    State__init_uniform(intron_end, "intron_end", 1, params->emission_table_4_UNIFORM);
    State__add_incoming(intron_end,            params->e1_e1,      intron_end);
    State__add_incoming(between_split_donor,   params->skip_do,    intron_end);
    State__add_incoming(match_donor_end,       params->do2_e1,      intron_end);
    former_intron = intron_end;

    struct Transition t;
    
    // where to start
    if (i == 0) {
      t.origin=first_intron->id;
      t.logodd=params->b1_b2;
      HMM__set_start(hmm, t);

      t.origin=match_acceptor->id;
      t.logodd=params->b1_acc;
      HMM__set_start(hmm, t);

      t.origin=between_acc_split->id;
      t.logodd=params->b1_bas;
      HMM__set_start(hmm, t);
    }

    // where to end
    if (i+1 == fasta->num_references) {
      t.origin = intron_end->id;
      t.logodd = params->e1_e2;
      HMM__set_end(hmm, t);

      t.origin = between_split_donor->id;
      t.logodd = params->bsd_e2;
      HMM__set_end(hmm, t);

      t.origin = match_donor_end->id;
      t.logodd = params->do2_e2;
      HMM__set_end(hmm, t);
    }

    if(params->dirty) {
      for(STATE_ID_T i=0; i < hmm->num_states; i++) {
        HMM__normalize(hmm, &hmm->states[i]);
      }
    }
  } // end for

#if !NONORMALIZE
  for (uint8_t i=1; i+1 < num_split_codons; i+=2) {
    State__add_incoming(split_codons[i],       params->intron_del,     split_codons[i+1]);
    HMM__normalize(hmm, split_codons[i]);
  }
#endif
  free(split_codons);

  return hmm;
}
