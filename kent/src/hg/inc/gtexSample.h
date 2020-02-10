/* gtexSample.h was originally generated by the autoSql program, which also 
 * generated gtexSample.c and gtexSample.sql.  This header links the database and
 * the RAM representation of objects. */

#ifndef GTEXSAMPLE_H
#define GTEXSAMPLE_H

#define GTEXSAMPLE_NUM_COLS 11

extern char *gtexSampleCommaSepFieldNames;

struct gtexSample
/* GTEx sample description */
    {
    struct gtexSample *next;  /* Next in singly linked list. */
    char *sampleId;	/* GTEx sample identifier */
    char *tissue;	/* Tissue name. Links to tissue table */
    char *donor;	/* GTEx subject identifier. Links to donor table */
    int autolysisScore;	/* Level of tissue self-digestion (0-3; none,mild,moderate,severe, -1 if unknown) */
    char *ischemicTime;	/* Time from tissue removal to preservation, in 4hr intervals */
    float rin;	/* RNA Integrity Number */
    char *collectionSites;	/* GTEx Biospecimen Source Site list */
    char *batchId;	/* Nucleic acid isolation batch ID */
    char *isolationType;	/* Type of nucleic acid isolation */
    char *isolationDate;	/* Date of nucleic acid isolation */
    char *pathNotes;	/* Pathology report notes */
    };

void gtexSampleStaticLoad(char **row, struct gtexSample *ret);
/* Load a row from gtexSample table into ret.  The contents of ret will
 * be replaced at the next call to this function. */

struct gtexSample *gtexSampleLoad(char **row);
/* Load a gtexSample from row fetched with select * from gtexSample
 * from database.  Dispose of this with gtexSampleFree(). */

struct gtexSample *gtexSampleLoadAll(char *fileName);
/* Load all gtexSample from whitespace-separated file.
 * Dispose of this with gtexSampleFreeList(). */

struct gtexSample *gtexSampleLoadAllByChar(char *fileName, char chopper);
/* Load all gtexSample from chopper separated file.
 * Dispose of this with gtexSampleFreeList(). */

#define gtexSampleLoadAllByTab(a) gtexSampleLoadAllByChar(a, '\t');
/* Load all gtexSample from tab separated file.
 * Dispose of this with gtexSampleFreeList(). */

struct gtexSample *gtexSampleCommaIn(char **pS, struct gtexSample *ret);
/* Create a gtexSample out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new gtexSample */

void gtexSampleFree(struct gtexSample **pEl);
/* Free a single dynamically allocated gtexSample such as created
 * with gtexSampleLoad(). */

void gtexSampleFreeList(struct gtexSample **pList);
/* Free a list of dynamically allocated gtexSample's */

void gtexSampleOutput(struct gtexSample *el, FILE *f, char sep, char lastSep);
/* Print out gtexSample.  Separate fields with sep. Follow last field with lastSep. */

#define gtexSampleTabOut(el,f) gtexSampleOutput(el,f,'\t','\n');
/* Print out gtexSample as a line in a tab-separated file. */

#define gtexSampleCommaOut(el,f) gtexSampleOutput(el,f,',',',');
/* Print out gtexSample as a comma separated list including final comma. */

/* -------------------------------- End autoSql Generated Code -------------------------------- */

#endif /* GTEXSAMPLE_H */

void gtexSampleCreateTable(struct sqlConnection *conn, char *table);
/* Create GTEx sample table of given name. */
