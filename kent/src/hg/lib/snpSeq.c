/* snpSeq.c was originally generated by the autoSql program, which also 
 * generated snpSeq.h and snpSeq.sql.  This module links the database and
 * the RAM representation of objects. */

/* Copyright (C) 2014 The Regents of the University of California 
 * See README in this or parent directory for licensing information. */

#include "common.h"
#include "linefile.h"
#include "dystring.h"
#include "jksql.h"
#include "snpSeq.h"


void snpSeqStaticLoad(char **row, struct snpSeq *ret)
/* Load a row from snpSeq table into ret.  The contents of ret will
 * be replaced at the next call to this function. */
{

ret->name = row[0];
ret->file_offset = sqlSigned(row[1]);
}

struct snpSeq *snpSeqLoad(char **row)
/* Load a snpSeq from row fetched with select * from snpSeq
 * from database.  Dispose of this with snpSeqFree(). */
{
struct snpSeq *ret;

AllocVar(ret);
ret->name = cloneString(row[0]);
ret->file_offset = sqlSigned(row[1]);
return ret;
}

struct snpSeq *snpSeqLoadAll(char *fileName) 
/* Load all snpSeq from a whitespace-separated file.
 * Dispose of this with snpSeqFreeList(). */
{
struct snpSeq *list = NULL, *el;
struct lineFile *lf = lineFileOpen(fileName, TRUE);
char *row[2];

while (lineFileRow(lf, row))
    {
    el = snpSeqLoad(row);
    slAddHead(&list, el);
    }
lineFileClose(&lf);
slReverse(&list);
return list;
}

struct snpSeq *snpSeqLoadAllByChar(char *fileName, char chopper) 
/* Load all snpSeq from a chopper separated file.
 * Dispose of this with snpSeqFreeList(). */
{
struct snpSeq *list = NULL, *el;
struct lineFile *lf = lineFileOpen(fileName, TRUE);
char *row[2];

while (lineFileNextCharRow(lf, chopper, row, ArraySize(row)))
    {
    el = snpSeqLoad(row);
    slAddHead(&list, el);
    }
lineFileClose(&lf);
slReverse(&list);
return list;
}

struct snpSeq *snpSeqCommaIn(char **pS, struct snpSeq *ret)
/* Create a snpSeq out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new snpSeq */
{
char *s = *pS;

if (ret == NULL)
    AllocVar(ret);
ret->name = sqlStringComma(&s);
ret->file_offset = sqlSignedComma(&s);
*pS = s;
return ret;
}

void snpSeqFree(struct snpSeq **pEl)
/* Free a single dynamically allocated snpSeq such as created
 * with snpSeqLoad(). */
{
struct snpSeq *el;

if ((el = *pEl) == NULL) return;
freeMem(el->name);
freez(pEl);
}

void snpSeqFreeList(struct snpSeq **pList)
/* Free a list of dynamically allocated snpSeq's */
{
struct snpSeq *el, *next;

for (el = *pList; el != NULL; el = next)
    {
    next = el->next;
    snpSeqFree(&el);
    }
*pList = NULL;
}

void snpSeqOutput(struct snpSeq *el, FILE *f, char sep, char lastSep) 
/* Print out snpSeq.  Separate fields with sep. Follow last field with lastSep. */
{
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->name);
if (sep == ',') fputc('"',f);
fputc(sep,f);
fprintf(f, "%d", el->file_offset);
fputc(lastSep,f);
}

/* -------------------------------- End autoSql Generated Code -------------------------------- */
