/* stsMapMouseNew.c was originally generated by the autoSql program, which also 
 * generated stsMapMouseNew.h and stsMapMouseNew.sql.  This module links the database and
 * the RAM representation of objects. */

/* Copyright (C) 2014 The Regents of the University of California 
 * See README in this or parent directory for licensing information. */

#include "common.h"
#include "linefile.h"
#include "dystring.h"
#include "jksql.h"
#include "stsMapMouseNew.h"


void stsMapMouseNewStaticLoad(char **row, struct stsMapMouseNew *ret)
/* Load a row from stsMapMouseNew table into ret.  The contents of ret will
 * be replaced at the next call to this function. */
{

ret->chrom = row[0];
ret->chromStart = sqlSigned(row[1]);
ret->chromEnd = sqlUnsigned(row[2]);
ret->name = row[3];
ret->score = sqlUnsigned(row[4]);
ret->identNo = sqlUnsigned(row[5]);
ret->ctgAcc = row[6];
ret->otherAcc = row[7];
ret->rhChrom = row[8];
ret->rhPos = atof(row[9]);
ret->rhLod = atof(row[10]);
ret->wigChr = row[11];
ret->wigPos = atof(row[12]);
ret->mgiChrom = row[13];
ret->mgiPos = atof(row[14]);
}

struct stsMapMouseNew *stsMapMouseNewLoad(char **row)
/* Load a stsMapMouseNew from row fetched with select * from stsMapMouseNew
 * from database.  Dispose of this with stsMapMouseNewFree(). */
{
struct stsMapMouseNew *ret;

AllocVar(ret);
ret->chrom = cloneString(row[0]);
ret->chromStart = sqlSigned(row[1]);
ret->chromEnd = sqlUnsigned(row[2]);
ret->name = cloneString(row[3]);
ret->score = sqlUnsigned(row[4]);
ret->identNo = sqlUnsigned(row[5]);
ret->ctgAcc = cloneString(row[6]);
ret->otherAcc = cloneString(row[7]);
ret->rhChrom = cloneString(row[8]);
ret->rhPos = atof(row[9]);
ret->rhLod = atof(row[10]);
ret->wigChr = cloneString(row[11]);
ret->wigPos = atof(row[12]);
ret->mgiChrom = cloneString(row[13]);
ret->mgiPos = atof(row[14]);
return ret;
}

struct stsMapMouseNew *stsMapMouseNewLoadAll(char *fileName) 
/* Load all stsMapMouseNew from a whitespace-separated file.
 * Dispose of this with stsMapMouseNewFreeList(). */
{
struct stsMapMouseNew *list = NULL, *el;
struct lineFile *lf = lineFileOpen(fileName, TRUE);
char *row[15];

while (lineFileRow(lf, row))
    {
    el = stsMapMouseNewLoad(row);
    slAddHead(&list, el);
    }
lineFileClose(&lf);
slReverse(&list);
return list;
}

struct stsMapMouseNew *stsMapMouseNewLoadAllByChar(char *fileName, char chopper) 
/* Load all stsMapMouseNew from a chopper separated file.
 * Dispose of this with stsMapMouseNewFreeList(). */
{
struct stsMapMouseNew *list = NULL, *el;
struct lineFile *lf = lineFileOpen(fileName, TRUE);
char *row[15];

while (lineFileNextCharRow(lf, chopper, row, ArraySize(row)))
    {
    el = stsMapMouseNewLoad(row);
    slAddHead(&list, el);
    }
lineFileClose(&lf);
slReverse(&list);
return list;
}

struct stsMapMouseNew *stsMapMouseNewCommaIn(char **pS, struct stsMapMouseNew *ret)
/* Create a stsMapMouseNew out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new stsMapMouseNew */
{
char *s = *pS;

if (ret == NULL)
    AllocVar(ret);
ret->chrom = sqlStringComma(&s);
ret->chromStart = sqlSignedComma(&s);
ret->chromEnd = sqlUnsignedComma(&s);
ret->name = sqlStringComma(&s);
ret->score = sqlUnsignedComma(&s);
ret->identNo = sqlUnsignedComma(&s);
ret->ctgAcc = sqlStringComma(&s);
ret->otherAcc = sqlStringComma(&s);
ret->rhChrom = sqlStringComma(&s);
ret->rhPos = sqlFloatComma(&s);
ret->rhLod = sqlFloatComma(&s);
ret->wigChr = sqlStringComma(&s);
ret->wigPos = sqlFloatComma(&s);
ret->mgiChrom = sqlStringComma(&s);
ret->mgiPos = sqlFloatComma(&s);
*pS = s;
return ret;
}

void stsMapMouseNewFree(struct stsMapMouseNew **pEl)
/* Free a single dynamically allocated stsMapMouseNew such as created
 * with stsMapMouseNewLoad(). */
{
struct stsMapMouseNew *el;

if ((el = *pEl) == NULL) return;
freeMem(el->chrom);
freeMem(el->name);
freeMem(el->ctgAcc);
freeMem(el->otherAcc);
freeMem(el->rhChrom);
freeMem(el->wigChr);
freeMem(el->mgiChrom);
freez(pEl);
}

void stsMapMouseNewFreeList(struct stsMapMouseNew **pList)
/* Free a list of dynamically allocated stsMapMouseNew's */
{
struct stsMapMouseNew *el, *next;

for (el = *pList; el != NULL; el = next)
    {
    next = el->next;
    stsMapMouseNewFree(&el);
    }
*pList = NULL;
}

void stsMapMouseNewOutput(struct stsMapMouseNew *el, FILE *f, char sep, char lastSep) 
/* Print out stsMapMouseNew.  Separate fields with sep. Follow last field with lastSep. */
{
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->chrom);
if (sep == ',') fputc('"',f);
fputc(sep,f);
fprintf(f, "%d", el->chromStart);
fputc(sep,f);
fprintf(f, "%u", el->chromEnd);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->name);
if (sep == ',') fputc('"',f);
fputc(sep,f);
fprintf(f, "%u", el->score);
fputc(sep,f);
fprintf(f, "%u", el->identNo);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->ctgAcc);
if (sep == ',') fputc('"',f);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->otherAcc);
if (sep == ',') fputc('"',f);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->rhChrom);
if (sep == ',') fputc('"',f);
fputc(sep,f);
fprintf(f, "%f", el->rhPos);
fputc(sep,f);
fprintf(f, "%f", el->rhLod);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->wigChr);
if (sep == ',') fputc('"',f);
fputc(sep,f);
fprintf(f, "%f", el->wigPos);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->mgiChrom);
if (sep == ',') fputc('"',f);
fputc(sep,f);
fprintf(f, "%f", el->mgiPos);
fputc(lastSep,f);
}

/* -------------------------------- End autoSql Generated Code -------------------------------- */

