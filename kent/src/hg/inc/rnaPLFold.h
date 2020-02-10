/* rnaPLFold.h was originally generated by the autoSql program, which also 
 * generated rnaPLFold.c and rnaPLFold.sql.  This header links the database and
 * the RAM representation of objects. */

/* Copyright (C) 2007 The Regents of the University of California 
 * See README in this or parent directory for licensing information. */

#ifndef RNAPLFOLD_H
#define RNAPLFOLD_H

#define RNAPLFOLD_NUM_COLS 5

struct rnaPLFold
/* RNA PL Fold -- RNA Partition Function Local Fold */
    {
    struct rnaPLFold *next;  /* Next in singly linked list. */
    char *chrom;	/* Reference sequence chromosome or scaffold */
    unsigned chromStart;	/* chromStart for reference marker */
    unsigned chromEnd;	/* chromEnd for last marker in list */
    unsigned span;	/* Number of positions covered by colorIndex values */
    char *colorIndex;	/* PLFold values encoded as indexes in a color array */
    };

void rnaPLFoldStaticLoad(char **row, struct rnaPLFold *ret);
/* Load a row from rnaPLFold table into ret.  The contents of ret will
 * be replaced at the next call to this function. */

struct rnaPLFold *rnaPLFoldLoad(char **row);
/* Load a rnaPLFold from row fetched with select * from rnaPLFold
 * from database.  Dispose of this with rnaPLFoldFree(). */

struct rnaPLFold *rnaPLFoldLoadAll(char *fileName);
/* Load all rnaPLFold from whitespace-separated file.
 * Dispose of this with rnaPLFoldFreeList(). */

struct rnaPLFold *rnaPLFoldLoadAllByChar(char *fileName, char chopper);
/* Load all rnaPLFold from chopper separated file.
 * Dispose of this with rnaPLFoldFreeList(). */

#define rnaPLFoldLoadAllByTab(a) rnaPLFoldLoadAllByChar(a, '\t');
/* Load all rnaPLFold from tab separated file.
 * Dispose of this with rnaPLFoldFreeList(). */

struct rnaPLFold *rnaPLFoldCommaIn(char **pS, struct rnaPLFold *ret);
/* Create a rnaPLFold out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new rnaPLFold */

void rnaPLFoldFree(struct rnaPLFold **pEl);
/* Free a single dynamically allocated rnaPLFold such as created
 * with rnaPLFoldLoad(). */

void rnaPLFoldFreeList(struct rnaPLFold **pList);
/* Free a list of dynamically allocated rnaPLFold's */

void rnaPLFoldOutput(struct rnaPLFold *el, FILE *f, char sep, char lastSep);
/* Print out rnaPLFold.  Separate fields with sep. Follow last field with lastSep. */

#define rnaPLFoldTabOut(el,f) rnaPLFoldOutput(el,f,'\t','\n');
/* Print out rnaPLFold as a line in a tab-separated file. */

#define rnaPLFoldCommaOut(el,f) rnaPLFoldOutput(el,f,',',',');
/* Print out rnaPLFold as a comma separated list including final comma. */

/* -------------------------------- End autoSql Generated Code -------------------------------- */

#define RNAPLFOLD_DIAGDARK "RNAplfold.diagDark"
#define RNAPLFOLD_DIAGDARK_DEF "50"
#define RNAPLFOLD_DIAGLIGHT "RNAplfold.diagLight"
#define RNAPLFOLD_DIAGLIGHT_DEF "10"
#define RNAPLFOLD_MAXBPDISTANCE "RNAplfold.maxBpDistance"
#define RNAPLFOLD_MAXBPDISTANCE_DEF "100"
#define RNAPLFOLD_INVERT "RNAplfold.invert"
#define RNAPLFOLD_INVERT_TOP "top"
#define RNAPLFOLD_INVERT_BUTTOM "buttom"
#define RNAPLFOLD_INVERT_DEF RNAPLFOLD_INVERT_BUTTOM

#endif /* RNAPLFOLD_H */

