# pgPhenoAssoc.sql was originally generated by the autoSql program, which also 
# generated pgPhenoAssoc.c and pgPhenoAssoc.h.  This creates the database representation of
# an object which can be loaded and saved from RAM in a fairly 
# automatic way.

#phenotypes from various databases for pgSnp tracks
CREATE TABLE pgPhenoAssoc (
    chrom varchar(255) not null,	# Chromosome
    chromStart int unsigned not null,	# Start position in chrom
    chromEnd int unsigned not null,	# End position in chrom
    name varchar(255) not null,	# Phenotype
    srcUrl longblob not null,	# link back to source database
              #Indices
    INDEX(chrom(12), chromStart)
);
