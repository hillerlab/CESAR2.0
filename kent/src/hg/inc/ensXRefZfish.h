/* ensXRefZfish.h was originally generated by the autoSql program, which also 
 * generated ensXRefZfish.c and ensXRefZfish.sql.  This header links the database and
 * the RAM representation of objects. */

/* Copyright (C) 2006 The Regents of the University of California 
 * See README in this or parent directory for licensing information. */

#ifndef ENSXREFZFISH_H
#define ENSXREFZFISH_H

#define ENSXREFZFISH_NUM_COLS 9

struct ensXRefZfish
/* Link from an Ensembl Transcript ID to other database IDs and description. */
    {
    struct ensXRefZfish *next;  /* Next in singly linked list. */
    char *ensGeneId;	/* Ensembl Transcript ID */
    char *zfinId;	/* ZFIN ID */
    char *uniProtId;	/* Unified UniProt protein accession */
    char *spDisplayId;	/* UniProt Display ID */
    char *geneId;	/* ZFIN Gene Symbol (formerly LocusLink) ID */
    char *geneSymbol;	/* Official ZFIN Gene Symbol */
    char *refSeq;	/* RefSeq DNA Accession */
    char *protAcc;	/* RefSeq Protein Accession */
    char *description;	/* Description */
    };

void ensXRefZfishStaticLoad(char **row, struct ensXRefZfish *ret);
/* Load a row from ensXRefZfish table into ret.  The contents of ret will
 * be replaced at the next call to this function. */

struct ensXRefZfish *ensXRefZfishLoad(char **row);
/* Load a ensXRefZfish from row fetched with select * from ensXRefZfish
 * from database.  Dispose of this with ensXRefZfishFree(). */

struct ensXRefZfish *ensXRefZfishLoadAll(char *fileName);
/* Load all ensXRefZfish from whitespace-separated file.
 * Dispose of this with ensXRefZfishFreeList(). */

struct ensXRefZfish *ensXRefZfishLoadAllByChar(char *fileName, char chopper);
/* Load all ensXRefZfish from chopper separated file.
 * Dispose of this with ensXRefZfishFreeList(). */

#define ensXRefZfishLoadAllByTab(a) ensXRefZfishLoadAllByChar(a, '\t');
/* Load all ensXRefZfish from tab separated file.
 * Dispose of this with ensXRefZfishFreeList(). */

struct ensXRefZfish *ensXRefZfishCommaIn(char **pS, struct ensXRefZfish *ret);
/* Create a ensXRefZfish out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new ensXRefZfish */

void ensXRefZfishFree(struct ensXRefZfish **pEl);
/* Free a single dynamically allocated ensXRefZfish such as created
 * with ensXRefZfishLoad(). */

void ensXRefZfishFreeList(struct ensXRefZfish **pList);
/* Free a list of dynamically allocated ensXRefZfish's */

void ensXRefZfishOutput(struct ensXRefZfish *el, FILE *f, char sep, char lastSep);
/* Print out ensXRefZfish.  Separate fields with sep. Follow last field with lastSep. */

#define ensXRefZfishTabOut(el,f) ensXRefZfishOutput(el,f,'\t','\n');
/* Print out ensXRefZfish as a line in a tab-separated file. */

#define ensXRefZfishCommaOut(el,f) ensXRefZfishOutput(el,f,',',',');
/* Print out ensXRefZfish as a comma separated list including final comma. */

/* -------------------------------- End autoSql Generated Code -------------------------------- */

#endif /* ENSXREFZFISH_H */

