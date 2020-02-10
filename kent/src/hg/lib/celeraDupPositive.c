/* celeraDupPositive.c was originally generated by the autoSql program, which also
 * generated celeraDupPositive.h and celeraDupPositive.sql.  This module links the database and the RAM
 * representation of objects. */

/* Copyright (C) 2014 The Regents of the University of California 
 * See README in this or parent directory for licensing information. */

#include "common.h"
#include "jksql.h"
#include "celeraDupPositive.h"


void celeraDupPositiveStaticLoad(char **row, struct celeraDupPositive *ret)
/* Load a row from celeraDupPositive table into ret.  The contents of ret will
 * be replaced at the next call to this function. */
{

ret->chrom = row[0];
ret->chromStart = sqlUnsigned(row[1]);
ret->chromEnd = sqlUnsigned(row[2]);
ret->name = row[3];
ret->fullName = row[4];
ret->fracMatch = atof(row[5]);
ret->bpAlign = atof(row[6]);
}

struct celeraDupPositive *celeraDupPositiveLoad(char **row)
/* Load a celeraDupPositive from row fetched with select * from celeraDupPositive
 * from database.  Dispose of this with celeraDupPositiveFree(). */
{
struct celeraDupPositive *ret;

AllocVar(ret);
ret->chrom = cloneString(row[0]);
ret->chromStart = sqlUnsigned(row[1]);
ret->chromEnd = sqlUnsigned(row[2]);
ret->name = cloneString(row[3]);
ret->fullName = cloneString(row[4]);
ret->fracMatch = atof(row[5]);
ret->bpAlign = atof(row[6]);
return ret;
}

struct celeraDupPositive *celeraDupPositiveCommaIn(char **pS, struct celeraDupPositive *ret)
/* Create a celeraDupPositive out of a comma separated string.
 * This will fill in ret if non-null, otherwise will
 * return a new celeraDupPositive */
{
char *s = *pS;

if (ret == NULL)
    AllocVar(ret);
ret->chrom = sqlStringComma(&s);
ret->chromStart = sqlUnsignedComma(&s);
ret->chromEnd = sqlUnsignedComma(&s);
ret->name = sqlStringComma(&s);
ret->fullName = sqlStringComma(&s);
ret->fracMatch = sqlSignedComma(&s);
ret->bpAlign = sqlSignedComma(&s);
*pS = s;
return ret;
}

void celeraDupPositiveFree(struct celeraDupPositive **pEl)
/* Free a single dynamically allocated celeraDupPositive such as created
 * with celeraDupPositiveLoad(). */
{
struct celeraDupPositive *el;

if ((el = *pEl) == NULL) return;
freeMem(el->chrom);
freeMem(el->name);
freez(pEl);
}

void celeraDupPositiveFreeList(struct celeraDupPositive **pList)
/* Free a list of dynamically allocated celeraDupPositive's */
{
struct celeraDupPositive *el, *next;

for (el = *pList; el != NULL; el = next)
    {
    next = el->next;
    celeraDupPositiveFree(&el);
    }
*pList = NULL;
}

void celeraDupPositiveOutput(struct celeraDupPositive *el, FILE *f, char sep, char lastSep)
/* Print out celeraDupPositive.  Separate fields with sep. Follow last field with lastSep. */
{
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->chrom);
if (sep == ',') fputc('"',f);
fputc(sep,f);
fprintf(f, "%u", el->chromStart);
fputc(sep,f);
fprintf(f, "%u", el->chromEnd);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->name);
if (sep == ',') fputc('"',f);
fputc(sep,f);
fprintf(f, "%s", el->fullName);
if (sep == ',') fputc('"',f);
fputc(sep,f);
fprintf(f, "%f", el->fracMatch);
fputc(sep,f);
fprintf(f, "%f", el->bpAlign);
fputc(lastSep,f);
}

