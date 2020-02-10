/* wgEncodeGencodeTranscriptSource.h was originally generated by the autoSql program, which also 
 * generated wgEncodeGencodeTranscriptSource.c and wgEncodeGencodeTranscriptSource.sql.  This header links the database and
 * the RAM representation of objects. */

/* Copyright (C) 2011 The Regents of the University of California 
 * See README in this or parent directory for licensing information. */

#ifndef WGENCODEGENCODETRANSCRIPTSOURCE_H
#define WGENCODEGENCODETRANSCRIPTSOURCE_H

#define WGENCODEGENCODETRANSCRIPTSOURCE_NUM_COLS 2

struct wgEncodeGencodeTranscriptSource
/* The source of Gencode transcript annotation */
    {
    struct wgEncodeGencodeTranscriptSource *next;  /* Next in singly linked list. */
    char *transcriptId;	/* GENCODE transcript identifier */
    char *source;	/* Source of transcript */
    };

void wgEncodeGencodeTranscriptSourceStaticLoad(char **row, struct wgEncodeGencodeTranscriptSource *ret);
/* Load a row from wgEncodeGencodeTranscriptSource table into ret.  The contents of ret will
 * be replaced at the next call to this function. */

struct wgEncodeGencodeTranscriptSource *wgEncodeGencodeTranscriptSourceLoad(char **row);
/* Load a wgEncodeGencodeTranscriptSource from row fetched with select * from wgEncodeGencodeTranscriptSource
 * from database.  Dispose of this with wgEncodeGencodeTranscriptSourceFree(). */

struct wgEncodeGencodeTranscriptSource *wgEncodeGencodeTranscriptSourceLoadAll(char *fileName);
/* Load all wgEncodeGencodeTranscriptSource from whitespace-separated file.
 * Dispose of this with wgEncodeGencodeTranscriptSourceFreeList(). */

struct wgEncodeGencodeTranscriptSource *wgEncodeGencodeTranscriptSourceLoadAllByChar(char *fileName, char chopper);
/* Load all wgEncodeGencodeTranscriptSource from chopper separated file.
 * Dispose of this with wgEncodeGencodeTranscriptSourceFreeList(). */

#define wgEncodeGencodeTranscriptSourceLoadAllByTab(a) wgEncodeGencodeTranscriptSourceLoadAllByChar(a, '\t');
/* Load all wgEncodeGencodeTranscriptSource from tab separated file.
 * Dispose of this with wgEncodeGencodeTranscriptSourceFreeList(). */

struct wgEncodeGencodeTranscriptSource *wgEncodeGencodeTranscriptSourceCommaIn(char **pS, struct wgEncodeGencodeTranscriptSource *ret);
/* Create a wgEncodeGencodeTranscriptSource out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new wgEncodeGencodeTranscriptSource */

void wgEncodeGencodeTranscriptSourceFree(struct wgEncodeGencodeTranscriptSource **pEl);
/* Free a single dynamically allocated wgEncodeGencodeTranscriptSource such as created
 * with wgEncodeGencodeTranscriptSourceLoad(). */

void wgEncodeGencodeTranscriptSourceFreeList(struct wgEncodeGencodeTranscriptSource **pList);
/* Free a list of dynamically allocated wgEncodeGencodeTranscriptSource's */

void wgEncodeGencodeTranscriptSourceOutput(struct wgEncodeGencodeTranscriptSource *el, FILE *f, char sep, char lastSep);
/* Print out wgEncodeGencodeTranscriptSource.  Separate fields with sep. Follow last field with lastSep. */

#define wgEncodeGencodeTranscriptSourceTabOut(el,f) wgEncodeGencodeTranscriptSourceOutput(el,f,'\t','\n');
/* Print out wgEncodeGencodeTranscriptSource as a line in a tab-separated file. */

#define wgEncodeGencodeTranscriptSourceCommaOut(el,f) wgEncodeGencodeTranscriptSourceOutput(el,f,',',',');
/* Print out wgEncodeGencodeTranscriptSource as a comma separated list including final comma. */

/* -------------------------------- End autoSql Generated Code -------------------------------- */

#endif /* WGENCODEGENCODETRANSCRIPTSOURCE_H */

