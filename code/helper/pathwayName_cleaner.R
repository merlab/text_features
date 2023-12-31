cleanNames <- function(v) {
  v <- gsub("_", " ", v)
  v <- gsub("[0-9]hr", "", v)
  v <- gsub("dependent", "dep.", v)
  v <- gsub("protein", "prot.", v)


  v <- gsub("NPC NPC", "NPC", v)

  v <- gsub("^REACTOME", "", v)
  v <- gsub("^WP", "", v)
  v <- gsub("^TCGA", "", v)
  v <- gsub("^MOOTHA", "", v)
  v <- gsub("^KEGG", "", v)
  v <- gsub("^NIKOLSKY", "", v)
  v <- gsub("^LINDGREN", "", v)
  v <- gsub("^BILD", "", v)
  v <- gsub("^GOBP", "", v)
  v <- gsub("^GOCC", "", v)
  v <- gsub("^GOMF", "", v)
  v <- gsub("^PID", "", v)
  v <- gsub("^QI", "", v)
  v <- gsub("^LEI", "", v)
  v <- gsub("^CHANDRAN", "", v)
  v <- gsub("MEDICUS", "", v)
  # v <- gsub("", "", v)
  # v <- gsub("", "", v)
  # v <- gsub("", "", v)
  # common abbreivation
  v <- gsub("EPITHELIAL TO MESENCHYMAL TRANSITION", "EMT", v)
  v <- gsub("NUCLEAR PORE COMPLEX", "NPC", v)
  v <- gsub("ELECTRON TRANSPORT ATP", "ET ATP", v)
  v <- gsub("ELECTRON TRANSPORT CHAIN", "ETC", v)
  v <- gsub("DOUBLE STRANDED BREAK", "DSB", v)
  v <- gsub("DOUBLE STRAND BREAK", "DSB", v)
  v <- gsub("CYTOCHROME P450", "CYP450", v)
  v <- gsub("GROWTH FACTOR RECEPTORS", "GFRs", v)
  v <- gsub("SECOND", "2nd", v)
  v <- gsub("G1 S", "G1/S", v)
  v <- gsub("G2 M", "G2/M", v)
  v <- gsub("RECEPTOR PROTEIN TYROSINE KINASE", "RTK", v)
  v <- gsub("RECEPTOR TYROSINE KINASE", "RTK", v)
  ###
  ###
  ###
  ###
  # v <- gsub("INVOLVED IN", "IN", v)
  # v <- gsub("THE ROLE OF", "ROLE OF", v)
  # v <- gsub(" THE ", " ", v)
  # v <- gsub(" AND ", " & ", v)
  # v <- gsub(" WITH ", " w/ ", v)
  ###
  ###
  ###
  ###
  # v <- gsub("PRODUCTION BY", "PROD.", v)
  # v <- gsub("PROTEINS", "PROT.", v)
  # v <- gsub("COUPLING", "CUPL.", v)
  # v <- gsub("PROTEIN", "PROT.", v)
  # v <- gsub(" & ", "+", v)
  # v <- gsub("", "", v)
  # v <- gsub("", "", v)
  # v <- gsub("", "", v)
  # v <- gsub("", "", v)
  # # shortening of documents
  # v <- gsub("MECHANISM", "MECH.", v)
  # v <- gsub("REGULATION OF", "REG.", v)
  # v <- gsub("CELLULAR", "CELL", v)
  # v <- gsub("NEGATIVE", "NEG.", v)
  # v <- gsub("SYSTEM", "SYS.", v)
  # v <- gsub("UBIQUITINATION", "UB.", v)
  # v <- gsub("REMOVAL", "RM.", v)
  # v <- gsub("REMOVE", "RM.", v)
  # # v <- gsub("EXPRESSION", "EXP.", v)
  # v <- gsub("ACTIVITY", "ACTY.", v)
  # v <- gsub("ACTIVATION", "ACTY.", v)
  # v <- gsub("FORMATION", "FORMN.", v)
  # v <- gsub("REPLICATION", "RPL.", v)
  # v <- gsub("REPLICATIVE", "RPL.", v)
  # v <- gsub("DEGRADATION", "DEG.", v)
  # v <- gsub("INITIATION", "INIT.", v)
  # v <- gsub("MOLECULES", "MOLEC.", v)
  # v <- gsub("MOLECULAR", "MOLEC.", v)
  # v <- gsub("MOLECULE", "MOLEC.", v)
  # v <- gsub("ASSOCIATED", "ASSOC.", v)
  # v <- gsub("ASSOCIATION", "ASSOC.", v)
  # v <- gsub("RESPONSE", "RESP.", v)
  # v <- gsub("SIGNALING", "SIG.", v)
  # v <- gsub("DAMAGE", "DMG.", v)
  # v <- gsub("RECOGNITION", "RECOG.", v)
  # v <- gsub("REPROGRAMMING", "REPROG.", v)
  # v <- gsub("PROGRAMMING", "PROG.", v)
  # v <- gsub("INTERACTIONS", "INT.", v)
  # v <- gsub("INTERACTION", "INT.", v)
  # v <- gsub("APOPTOSISRELATED NETWORK", "APOPTOSIS RELATED NTW.", v)
  # v <- gsub("ALTERED", "ALT.", v)
  # v <- gsub("NETWORK", "NTW.", v)
  # v <- gsub("METABOLIC", "METAB.", v)
  # v <- gsub("METABOLISM", "METAB.", v)
  # v <- gsub("CHECKPOINTS", "CKPT.", v)
  # v <- gsub("CHECKPOINT", "CKPT.", v)
  # v <- gsub("MAINTENANCE", "MAINT.", v)
  # v <- gsub("SYNTHESIS", "SYN.", v)
  # # v <- gsub("COMPLEX", "CMPLX.", v)
  # v <- gsub("EXTENSION", "EXT.", v)
  # v <- gsub("PROGRESSION", "PROG.", v)
  # v <- gsub("PHOSPHORYLATION", "PY.", v)
  # v <- gsub("STABILITY", "STB.", v)
  # v <- gsub("TARGETTED", "TARGET", v)
  # v <- gsub("TELOMERES", "telomere", v)
  # lexical simplication
  # v <- gsub(" PATHWAY ", " ", v)
  # v <- gsub(" THROUGH ", " via ", v)
  # v <- gsub(" VIA ", " via ", v)
  # v <- gsub(" MEDIATED ", " med. ", v)
  # v <- gsub(" OF ", " of ", v)
  # v <- gsub(" TO ", " to ", v)
  # v <- gsub(" IN ", " in ", v)
  # v <- gsub(" PRE ", " pre ", v)
  # v <- gsub(" DUE ", " due ", v)
  # v <- gsub(" due to ", " from ", v)
  # v <- gsub("  ", "  ", v)
  # v <- gsub("  ", "  ", v)
  # v <- gsub(" DERIVED FROM ", " from ", v)
  #
  # v <- gsub("OVARIAN", "OV.", v)
  # # v <- gsub("", "", v)
  # # v <- gsub("", "", v)
  # # v <- gsub("", "", v)
  # v <- gsub("TRANSITION", "TRANS.", v)
  # v <- gsub("TRANSDUCTION", "TRANS.", v)
  # v <- gsub("MITOTIC G1/S", "G1/S", v)
  v <- gsub("^ ", "", v)
  #
  return(v)
}
rmExtra <- function(mat) {
  mat <- mat[!grepl("infection", mat$pathway, ignore.case = TRUE), ]
  mat <- mat[!grepl("alzheimers", mat$pathway, ignore.case = TRUE), ]
  mat <- mat[!grepl("biocarta", mat$pathway, ignore.case = TRUE), ]
  mat <- mat[!grepl("infect", mat$pathway, ignore.case = TRUE), ]
  mat <- mat[!grepl("covid", mat$pathway, ignore.case = TRUE), ]
  mat <- mat[!grepl("sars", mat$pathway, ignore.case = TRUE), ]
  mat <- mat[!grepl("influencsa", mat$pathway, ignore.case = TRUE), ]
  mat <- mat[!grepl("parkinson", mat$pathway, ignore.case = TRUE), ]
  mat <- mat[!grepl("huntingtons", mat$pathway, ignore.case = TRUE), ]
  mat <- mat[!grepl("CARBOHYDRATE", mat$pathway, ignore.case = TRUE), ]
  mat <- mat[!grepl("kidney", mat$pathway, ignore.case = TRUE), ]
  mat <- mat[!grepl("heat", mat$pathway, ignore.case = TRUE), ]
  mat <- mat[!grepl("copy", mat$pathway, ignore.case = TRUE), ]
  mat <- mat[!grepl("^metastasis", mat$pathway, ignore.case = TRUE), ]
  mat <- mat[!grepl("^innate", mat$pathway, ignore.case = TRUE), ]
  mat <- mat[!grepl("POSITIVE REGULATION OF MOLECULAR FUNCTION", mat$pathway, ignore.case = TRUE), ]
  # mat <- mat[!grepl("", mat$pathway, ignore.case = TRUE), ]
  # mat <- mat[!grepl("lung", mat$pathway, ignore.case = TRUE), ]
  return(mat)
}
