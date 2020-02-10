table fishClones
"Describes the positions of fishClones in the assembly"
   (
   string chrom;                   "Reference sequence chromosome or scaffold"
   uint   chromStart;              "Start position in chromosome"
   uint   chromEnd;                "End position in chromosome"
   string name;                    "Name of clone"
   uint score;                     "Score from 0-1000"
   uint placeCount;                "Number of FISH clone mappings. Some mappings may be from non-sequence based methods"
   string[placeCount] bandStarts;  "Starting band of FISH clone mapping"
   string[placeCount] bandEnds;    "Ending band of FISH clone mapping"
   string[placeCount] labs;        "Lab where clone FISH'd"
   string placeType;               "How clone was placed on the sequence assembly"
   uint accCount;                  "Number of accessions associated with the clone"
   string[accCount] accNames;      "Accession associated with clone"
   uint stsCount;                  "Number of STS markers associated with this clone"
   string[stsCount] stsNames;      "Names of STS  markers"
   uint beCount;                   "Number of BAC end sequences associated with this clone"
   string[beCount] beNames;        "Accessions of BAC ends"
   )