# Generates the SRL_<LLM> datasets: a random 300-row sample of the MSLQ
# responses each "valid" (structured) LLM produced in the Computers in Human
# Behavior (2025) study. The two incoherent generators (ChatGPT, LeChat, whose
# items carry no estimable partial-correlation structure) are excluded.
# Source CSV is local-only; this script documents how the shipped data was made.
src <- "/Users/mohammedsaqr/Downloads/rstudio-export (18) 2/mohammed_all.csv"
d <- read.csv(src, check.names = FALSE)

# Ship the 44 MSLQ items only. The construct composites in the source CSV are
# exactly the per-construct item means (redundant and perfectly collinear), so
# they are dropped; recompute them from the item suffixes when needed.
items <- grep("^Q_", names(d), value = TRUE)        # 44 MSLQ Likert items
valid <- c("GPT", "Gemini", "Claude", "Mistral", "LLaMa")

set.seed(2025)
for (llm in valid) {
  rows <- which(d$Source == llm)
  take <- sort(sample(rows, 300L))
  df <- d[take, items]
  rownames(df) <- NULL
  assign(paste0("SRL_", llm), df)
  save(list = paste0("SRL_", llm),
       file = file.path("data", paste0("SRL_", llm, ".rda")), compress = "xz")
}
