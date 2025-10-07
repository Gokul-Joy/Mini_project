
##########################
# Step 1: Literature & Database Cross-check
##########################

# Install/load required packages
packages <- c("easyPubMed", "jsonlite", "httr", "dplyr", "ggplot2", "stringr")
invisible(lapply(packages, function(pkg) {
  if (!require(pkg, character.only = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}))

# ---- INPUT ----
biomarkers <- c("ARAP1"  , "CMTM8"  , "CXorf66" ,"FAM83A" , "FAM83F"  ,"KLF1"   , "KRT19" ,
                "LRRC1"   ,"MYEOV"  , "OR51F2" ,
                "PCDH1"   ,"PERP"  ,  "S100A5" , "SPIRE2"  ,"USH1C" ) 

disease_keyword <- "pancreatic cancer"    


# ---- 1. PubMed search counts ----
get_pubmed_count <- function(gene, disease) {
  query <- paste0(gene, " AND ", disease, "[Title/Abstract]")
  res <- easyPubMed::get_pubmed_ids(query)
  data.frame(
    gene = gene,
    disease = disease,
    pub_count = res$Count
  )
}

pub_counts <- do.call(rbind, lapply(biomarkers, get_pubmed_count, disease_keyword))
print(pub_counts)





#PubMed yearly trends===========================================================
library(easyPubMed)
library(dplyr)
library(ggplot2)

# Custom yearly trend function for any easyPubMed version
get_pubmed_trend_safe <- function(gene, disease) {
  tryCatch({
    query <- paste0(gene, " AND ", disease, "[Title/Abstract]")
    res <- get_pubmed_ids(query)
    
    # If no results
    if (is.null(res$Count) || as.numeric(res$Count) == 0) {
      message("No PubMed results for: ", gene)
      return(tibble(year = integer(), n = integer(), gene = character()))
    }
    
    # Fetch articles & convert to data frame
    df <- fetch_pubmed_data(res, format = "xml") %>%
      article_to_df(max_chars = 500)
    
    # Extract and clean years
    df$year <- as.numeric(substr(df$year, 1, 4))
    
    # Summarise by year
    trend <- df %>%
      filter(!is.na(year)) %>%
      group_by(year) %>%
      summarise(n = n(), .groups = "drop") %>%
      mutate(gene = gene)
    
    Sys.sleep(0.3) # avoid rate limiting
    return(trend)
  },
  error = function(e) {
    message("Error for ", gene, ": ", e$message)
    return(tibble(year = integer(), n = integer(), gene = character()))
  })
}



pub_trends <- bind_rows(
  lapply(biomarkers, get_pubmed_trend_safe, disease = disease_keyword)
)

# Plot trends
ggplot(pub_trends, aes(x = year, y = n, color = gene)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  theme_minimal(base_size = 14) +
  labs(
    title = paste("Publication Trends for Biomarkers in", disease_keyword),
    x = "Year",
    y = "Publications per Year",
    color = "Gene"
  )















# ---- 3. DisGeNET API query ----
# API docs: https://www.disgenet.org/api/
# (Requires free registration + API key)
disgenet_key <- "YOUR_DISGENET_API_KEY"

get_disgenet <- function(gene) {
  url <- paste0("https://www.disgenet.org/api/gda/gene/", gene, "?source=ALL")
  res <- GET(url, add_headers(Authorization = paste("Bearer", disgenet_key)))
  if (status_code(res) == 200) {
    df <- fromJSON(content(res, "text"))
    return(df[, c("gene_symbol", "disease_name", "score")])
  } else {
    return(data.frame())
  }
}

disgenet_results <- do.call(rbind, lapply(biomarkers, get_disgenet))
print(disgenet_results)

# ---- 4. GeneCards & OMIM notes ----
# GeneCards doesn't have an API, but you can store URLs for manual lookup:
genecards_urls <- paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=", biomarkers)
omim_urls <- paste0("https://www.omim.org/search?search=", biomarkers)

urls_df <- data.frame(
  gene = biomarkers,
  genecards = genecards_urls,
  omim = omim_urls
)
print(urls_df)

# ---- OUTPUT ----
# Save all results
write.csv(pub_counts, "pubmed_counts.csv", row.names = FALSE)
write.csv(pub_trends, "pubmed_trends.csv", row.names = FALSE)
write.csv(disgenet_results, "disgenet_results.csv", row.names = FALSE)
write.csv(urls_df, "gene_urls.csv", row.names = FALSE)
